/* Histology Measure: feature finder (canals, osteons, cavities).
 * Candidates are connected regions darker (or brighter) than the tissue around them, with contrast normalised locally
 * inside the tissue so the section edge does not light up. The cartoon beside the image shows them colour-coded and
 * moves with it. The measurer labels them by clicking, by key (1-6, X), or by "Find similar" from labelled examples.
 * Labelled features are saved as pixel runs in the export; a mask PNG beside the export holds the whole segmentation
 * (labelled and unlabelled), so an export can be reloaded and relabelled later.
 */
(function () {
"use strict";
const HM = window.HM; if (!HM) return;
const { S, $, h } = HM;

const CLASSES = [
  ["canal", "Vascular canal", "#dc2626"],
  ["primary_osteon", "Primary osteon", "#e08600"],
  ["secondary_osteon", "Secondary osteon", "#7c3aed"],
  ["resorption_cavity", "Resorption cavity", "#0e9fbf"],
  ["artefact", "Crack / artefact", "#475569"],
  ["other", "Other", "#db2777"],
];
const NONE = "none";                                   // "not a feature": kept, and used as a negative example
const NAME = Object.fromEntries(CLASSES.map(([k, l]) => [k, l])); NAME[NONE] = "Not a feature";
const HEX = Object.fromEntries(CLASSES.map(([k, , c]) => [k, c])); HEX[NONE] = "#d9d4c9";
const DEFAULTS = { light: null, polarity: "dark", window_px: 8, threshold: 1.0, min_area_px: 12, context_px: 30, match: 0.9 };
const THR = [3, 0.4];                                  // sensitivity slider: 0 -> threshold 3 SD, 100 -> 0.4 SD
const MASK_ENCODING = "value = R + 256*G + 65536*B; 0 = no feature. Values in objects[].mask_value are labelled; other values are unlabelled candidates.";
const rgb = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16));
const mix = (a, b, t) => a.map((v, i) => Math.round(v * (1 - t) + b[i] * t));
const C_OUT = rgb("#2a2a2e"), C_OFF = rgb("#6a6a72"), C_TIS = rgb("#f1ece1"), C_CAND = rgb("#8d929b");
const C_LAB = Object.fromEntries(Object.entries(HEX).map(([k, c]) => [k, rgb(c)]));
const C_SUG = Object.fromEntries(Object.entries(HEX).map(([k, c]) => [k, mix(rgb(c), C_TIS, 0.55)]));

// Working state for the open image. Only doc.features (settings + labelled objects) is saved.
const F = { url: null, px: null, tissue: null, tm: 128, ts: 40, score: null, scoreKey: "", region: null, regionKind: "",
            regionPx: 0, regionSig: "", lab: null, objs: [], suggest: new Map(), hover: -1, on: false,
            raster: null, rctx: null, rimg: null, busy: false, timer: null, fromMask: false };

function fs(d = S.doc) {
  if (!d.features) d.features = { settings: {}, objects: [] };
  const st = d.features.settings;
  for (const k in DEFAULTS) if (st[k] == null) st[k] = DEFAULTS[k];
  if (!st.light) st.light = /pol|xpl|cross|lambda/i.test(String((d.setup && d.setup.light) || d.task.light || "")) ? "xpl" : "ppl";
  return d.features;
}
const current = () => S.doc && S.img && F.px && F.url === S.img.url;

// ---------- pixels ----------
async function ensurePixels() {
  if (current()) return true;
  const url = S.img.url, im = new Image(); im.src = url;
  try { await im.decode(); } catch (e) { HM.status("Pixels unreadable"); return false; }
  if (!S.img || S.img.url !== url) return false;
  const w = im.naturalWidth, hh = im.naturalHeight, N = w * hh;
  const cv = document.createElement("canvas"); cv.width = w; cv.height = hh;
  const g = cv.getContext("2d", { willReadFrequently: true }); g.drawImage(im, 0, 0);
  const rgba = g.getImageData(0, 0, w, hh).data, L = new Float32Array(N);
  for (let i = 0, j = 0; i < N; i++, j += 4) L[i] = 0.299 * rgba[j] + 0.587 * rgba[j + 1] + 0.114 * rgba[j + 2];
  reset(); Object.assign(F, { url, px: { w, h: hh, N, rgba, L } });
  F.tissue = tissueMask();
  let s = 0, s2 = 0, n = 0;
  for (let i = 0; i < N; i += 7) if (F.tissue[i]) { s += L[i]; s2 += L[i] * L[i]; n++; }
  F.tm = n ? s / n : 128; F.ts = n ? Math.sqrt(Math.max(s2 / n - F.tm * F.tm, 1)) : 40;
  return true;
}
function reset() {
  Object.assign(F, { url: null, px: null, tissue: null, score: null, scoreKey: "", region: null, regionSig: "", lab: null,
                     objs: [], suggest: new Map(), hover: -1, raster: null, rctx: null, rimg: null, fromMask: false });
}
// Slide background: pixels close to the border colour and connected to the border (white in plain light, black under
// crossed polars). Dark canals inside the section are not connected to the border, so they stay tissue.
function tissueMask() {
  const { w, h: hh, N, rgba } = F.px, R = [], G = [], B = [];
  const take = (x, y) => { const j = 4 * (y * w + x); R.push(rgba[j]); G.push(rgba[j + 1]); B.push(rgba[j + 2]); };
  for (let x = 0; x < w; x += 2) { take(x, 0); take(x, hh - 1); }
  for (let y = 0; y < hh; y += 2) { take(0, y); take(w - 1, y); }
  const med = (a) => a.sort((p, q) => p - q)[a.length >> 1], bg = [med(R), med(G), med(B)];
  const out = new Uint8Array(N).fill(1), q = new Int32Array(N); let qn = 0;
  const seed = (i) => { if (!out[i]) return; const j = 4 * i;
    if (Math.abs(rgba[j] - bg[0]) + Math.abs(rgba[j + 1] - bg[1]) + Math.abs(rgba[j + 2] - bg[2]) < 45) { out[i] = 0; q[qn++] = i; } };
  for (let x = 0; x < w; x++) { seed(x); seed((hh - 1) * w + x); }
  for (let y = 0; y < hh; y++) { seed(y * w); seed(y * w + w - 1); }
  for (let k = 0; k < qn; k++) { const i = q[k], x = i % w;
    if (x > 0) seed(i - 1); if (x < w - 1) seed(i + 1); if (i >= w) seed(i - w); if (i < N - w) seed(i + w); }
  // Large enclosed slide-coloured regions (a medullary cavity joined to the edge only by hairline cracks) are not tissue.
  const near = new Uint8Array(N);
  for (let i = 0; i < N; i++) if (out[i]) { const j = 4 * i;
    near[i] = Math.abs(rgba[j] - bg[0]) + Math.abs(rgba[j + 1] - bg[1]) + Math.abs(rgba[j + 2] - bg[2]) < 45 ? 1 : 0; }
  for (const r of components(near, w, hh)) if (runArea(r) > 0.002 * N) for (let k = 0; k < r.length; k += 3) out.fill(0, r[k] * w + r[k + 1], r[k] * w + r[k + 2] + 1);
  return out;
}
// Separable box blur (mean over a (2r+1)^2 window, shrunk at the image edge).
function blur(src, w, hh, r) {
  const tmp = new Float32Array(w * hh), out = new Float32Array(w * hh), cs = new Float64Array(Math.max(w, hh) + 1);
  for (let y = 0; y < hh; y++) { const o = y * w;
    for (let x = 0; x < w; x++) cs[x + 1] = cs[x] + src[o + x];
    for (let x = 0; x < w; x++) { const a = Math.max(0, x - r), b = Math.min(w, x + r + 1); tmp[o + x] = (cs[b] - cs[a]) / (b - a); } }
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < hh; y++) cs[y + 1] = cs[y] + tmp[y * w + x];
    for (let y = 0; y < hh; y++) { const a = Math.max(0, y - r), b = Math.min(hh, y + r + 1); out[y * w + x] = (cs[b] - cs[a]) / (b - a); } }
  return out;
}
// Local contrast: how much darker (or brighter) each pixel is than the tissue within window_px, in local SDs.
function computeScore(st) {
  const key = st.polarity + "|" + st.window_px; if (F.score && F.scoreKey === key) return;
  const { w, h: hh, N, L } = F.px, m = F.tissue, Ls = blur(L, w, hh, 1);
  const mm = new Float32Array(N), m1 = new Float32Array(N), m2 = new Float32Array(N);
  for (let i = 0; i < N; i++) if (m[i]) { mm[i] = 1; m1[i] = Ls[i]; m2[i] = Ls[i] * Ls[i]; }
  const r = Math.max(2, Math.round(st.window_px)), bm = blur(mm, w, hh, r), b1 = blur(m1, w, hh, r), b2 = blur(m2, w, hh, r);
  const z = new Float32Array(N), pol = st.polarity;
  for (let i = 0; i < N; i++) {
    if (!m[i] || bm[i] < 0.05) { z[i] = -99; continue; }
    const mu = b1[i] / bm[i], v = b2[i] / bm[i] - mu * mu, d = (mu - Ls[i]) / (Math.sqrt(Math.max(v, 0)) + 4);
    z[i] = pol === "bright" ? -d : pol === "either" ? Math.abs(d) : d;
  }
  F.score = z; F.scoreKey = key;
}
// Region searched: traced Areas if any; else the cortex inside the periosteum (minus the medullary cavity) when both
// go most of the way round the centre; else all tissue. Ignore boxes are always left out.
function regionSig() { const d = S.doc;
  return JSON.stringify([d.centre, d.curves.filter((c) => ["count_area", "periosteum", "medullary", "ignore"].includes(c.kind)).map((c) => [c.kind, c.segments])]); }
const ring = (c, z) => c.segments.flatMap((s) => HM.resample(s, 4)).sort((a, b) => HM.angle(a, z) - HM.angle(b, z));
function computeRegion() {
  const d = S.doc, { w, h: hh, N } = F.px, z = d.centre;
  const areas = d.curves.filter((c) => c.kind === "count_area" && c.segments[0] && c.segments[0].length > 2);
  const peri = d.curves.find((c) => c.kind === "periosteum"), med = d.curves.find((c) => c.kind === "medullary");
  const round = (c) => c && z && ((HM.angular(c, z) || {}).coverage_deg || 0) >= 300;
  const ignore = d.curves.filter((c) => c.kind === "ignore" && c.segments[0] && c.segments[0].length > 2);
  const R = new Uint8Array(N); let kind = "tissue";
  const cv = document.createElement("canvas"); cv.width = w; cv.height = hh;
  const g = cv.getContext("2d", { willReadFrequently: true }); g.fillStyle = "#fff";
  const poly = (pts) => { g.beginPath(); pts.forEach((p, i) => (i ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1]))); g.closePath(); g.fill(); };
  if (areas.length) { kind = "areas"; areas.forEach((c) => poly(c.segments[0])); }
  else if (round(peri)) { kind = "cortex"; poly(ring(peri, z)); }
  else g.fillRect(0, 0, w, hh);
  g.globalCompositeOperation = "destination-out";
  if (kind === "cortex" && round(med)) poly(ring(med, z));
  ignore.forEach((c) => poly(c.segments[0]));
  const a = g.getImageData(0, 0, w, hh).data;
  for (let i = 0; i < N; i++) R[i] = a[4 * i + 3] > 127 && F.tissue[i] ? 1 : 0;
  let n = 0; for (let i = 0; i < N; i++) n += R[i];
  Object.assign(F, { region: R, regionKind: kind, regionPx: n });
}
function prepare(st) { computeScore(st); const sig = regionSig(); if (!F.region || F.regionSig !== sig) { computeRegion(); F.regionSig = sig; } }

// ---------- objects ----------
// Connected components (8-neighbour) by runs; each object is a flat list of runs [y, x0, x1, ...].
function components(mask, w, hh) {
  const ry = [], r0 = [], r1 = [], par = [];
  const find = (a) => { while (par[a] !== a) { par[a] = par[par[a]]; a = par[a]; } return a; };
  const unite = (a, b) => { a = find(a); b = find(b); if (a !== b) par[Math.max(a, b)] = Math.min(a, b); };
  let ps = 0, pe = 0;
  for (let y = 0; y < hh; y++) { const o = y * w, rs = ry.length; let p = ps;
    for (let x = 0; x < w;) { if (!mask[o + x]) { x++; continue; }
      const x0 = x; while (x < w && mask[o + x]) x++; const x1 = x - 1, id = ry.length;
      ry.push(y); r0.push(x0); r1.push(x1); par.push(id);
      while (p < pe && r1[p] < x0 - 1) p++;
      for (let q = p; q < pe && r0[q] <= x1 + 1; q++) unite(id, q); }
    ps = rs; pe = ry.length; }
  const groups = new Map();
  for (let i = 0; i < ry.length; i++) { const r = find(i); let g = groups.get(r); if (!g) groups.set(r, (g = [])); g.push(ry[i], r0[i], r1[i]); }
  return [...groups.values()];
}
// A connected network (reticular canals) larger than one tile is cut on a T-px grid, so each piece can be labelled.
function split(runs, T) {
  let x0 = 1e9, x1 = -1, y0 = 1e9, y1 = -1;
  for (let k = 0; k < runs.length; k += 3) { y0 = Math.min(y0, runs[k]); y1 = Math.max(y1, runs[k]); x0 = Math.min(x0, runs[k + 1]); x1 = Math.max(x1, runs[k + 2]); }
  if (x1 - x0 < T && y1 - y0 < T) return [runs];
  const cells = new Map();
  for (let k = 0; k < runs.length; k += 3) { const y = runs[k], cy = Math.floor(y / T);
    for (let a = runs[k + 1]; a <= runs[k + 2];) { const cx = Math.floor(a / T), b = Math.min(runs[k + 2], (cx + 1) * T - 1), key = cy * 100000 + cx;
      if (!cells.has(key)) cells.set(key, []); cells.get(key).push(y, a, b); a = b + 1; } }
  return [...cells.values()];
}
const runArea = (r) => { let n = 0; for (let k = 0; k < r.length; k += 3) n += r[k + 2] - r[k + 1] + 1; return n; };
function paintRuns(lab, r, v) { const w = F.px.w; for (let k = 0; k < r.length; k += 3) lab.fill(v, r[k] * w + r[k + 1], r[k] * w + r[k + 2] + 1); }
// Shape, tone, colour, orientation relative to the centre, and the tissue ring around the object (what separates an
// osteon's canal from a plain canal). Ring statistics are rotation-invariant.
function stats(r) {
  const { w, h: hh, rgba, L } = F.px, z = F.score, c0 = S.doc.centre, ctx = fs().settings.context_px;
  let n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sL = 0, sZ = 0, sr = 0, sg = 0, sb = 0, bx0 = 1e9, by0 = 1e9, bx1 = -1, by1 = -1;
  for (let k = 0; k < r.length; k += 3) { const y = r[k], a = r[k + 1], b = r[k + 2];
    if (y < by0) by0 = y; if (y > by1) by1 = y; if (a < bx0) bx0 = a; if (b > bx1) bx1 = b;
    for (let x = a; x <= b; x++) { const i = y * w + x, j = 4 * i; n++; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      sL += L[i]; sZ += Math.max(-10, z[i]); sr += rgba[j]; sg += rgba[j + 1]; sb += rgba[j + 2]; } }
  const cx = sx / n, cy = sy / n, vxx = sxx / n - cx * cx + 1 / 12, vyy = syy / n - cy * cy + 1 / 12, vxy = sxy / n - cx * cy;
  const tr = vxx + vyy, det = vxx * vyy - vxy * vxy, disc = Math.sqrt(Math.max(tr * tr / 4 - det, 0));
  const l1 = tr / 2 + disc, l2 = Math.max(tr / 2 - disc, 1e-6), theta = 0.5 * Math.atan2(2 * vxy, vxx - vyy);
  const ringM = [], ringS = [];
  for (let q = 1; q <= 4; q++) { const rad = q * ctx / 4; let s = 0, s2 = 0;
    for (let t = 0; t < 16; t++) { const a = t * Math.PI / 8;
      const x = Math.min(w - 1, Math.max(0, Math.round(cx + rad * Math.cos(a)))), y = Math.min(hh - 1, Math.max(0, Math.round(cy + rad * Math.sin(a))));
      const v = (L[y * w + x] - F.tm) / F.ts; s += v; s2 += v * v; }
    ringM.push(s / 16); ringS.push(Math.sqrt(Math.max(s2 / 16 - (s / 16) ** 2, 0))); }
  const sum = sr + sg + sb || 1;
  return { n, cx, cy, bbox: [bx0, by0, bx1, by1], elong: 1 - Math.sqrt(l2 / l1), theta,
           radial: c0 ? Math.abs(Math.cos(theta - Math.atan2(cy - c0[1], cx - c0[0]))) : null,
           meanL: sL / n, meanZ: sZ / n, chr: sr / sum, chg: sg / sum, fill: n / ((bx1 - bx0 + 1) * (by1 - by0 + 1)), ringM, ringS };
}
const vecOf = (s) => [Math.log(s.n), s.elong, s.radial == null ? 0.5 : s.radial, (s.meanL - F.tm) / F.ts, s.meanZ,
  s.chr * 10, s.chg * 10, s.fill, ...s.ringM, ...s.ringS];
function standardize(objs) {
  if (!objs.length) return; const V = objs.map((o) => vecOf(o.st)), D = V[0].length, mu = [], sd = [];
  for (let j = 0; j < D; j++) { let s = 0, s2 = 0; for (const v of V) { s += v[j]; s2 += v[j] * v[j]; }
    mu.push(s / V.length); sd.push(Math.sqrt(Math.max(s2 / V.length - (s / V.length) ** 2, 1e-6))); }
  objs.forEach((o, i) => { o.v = V[i].map((x, j) => (x - mu[j]) / sd[j]); });
}
const r1 = (v) => (v == null ? null : Math.round(v * 10) / 10), r3 = (v) => (v == null ? null : Math.round(v * 1000) / 1000);
const exportStats = (s) => ({ centroid: [r1(s.cx), r1(s.cy)], area_px: s.n, bbox: s.bbox, elongation: r3(s.elong),
  orientation_deg: r1(s.theta * 180 / Math.PI), radial_cos: r3(s.radial), mean_L: r1(s.meanL) });

function finish(lab, objs) {
  for (const o of objs) { o.st = stats(o.runs); if (o.ref) Object.assign(o.ref, exportStats(o.st)); }
  standardize(objs);
  Object.assign(F, { lab, objs, suggest: new Map(), hover: -1 });
  paintAll(); HM.renderAll(); renderPanel();
}
async function detect() {
  if (!S.doc || !S.img || F.busy) return;
  F.busy = true; renderPanel(); HM.status("Finding features…", true);
  await new Promise((r) => setTimeout(r, 30));
  try {
    if (!(await ensurePixels())) return;
    const st = fs().settings; prepare(st);
    const { w, h: hh, N } = F.px, lab = new Int32Array(N).fill(-1), objs = [];
    for (const o of fs().objects) { paintRuns(lab, o.runs, objs.length); objs.push({ ref: o, label: o.label, runs: o.runs }); }
    const mask = new Uint8Array(N), z = F.score, R = F.region, thr = st.threshold;
    for (let i = 0; i < N; i++) mask[i] = R[i] && lab[i] < 0 && z[i] > thr ? 1 : 0;
    const T = Math.max(24, 4 * st.window_px); let nets = 0;
    for (const runs of components(mask, w, hh)) {
      if (runArea(runs) < st.min_area_px) continue;
      const pieces = split(runs, T); if (pieces.length > 1) nets++;
      for (const r of pieces) { if (runArea(r) < Math.min(st.min_area_px, 4)) continue; paintRuns(lab, r, objs.length); objs.push({ ref: null, label: null, runs: r }); } }
    F.fromMask = false; finish(lab, objs);
    HM.status(`${objs.filter((o) => !o.ref).length} candidates` + (nets ? ` · ${nets} networks tiled` : ""));
  } finally { F.busy = false; renderPanel(); }
}
function schedule() { clearTimeout(F.timer); if (F.lab) F.timer = setTimeout(detect, 350); }

// ---------- labelling ----------
function setLabel(i, label, source = "manual") {
  const o = F.objs[i]; if (!o) return; const f = fs();
  if (label == null) { if (o.ref) f.objects.splice(f.objects.indexOf(o.ref), 1); o.ref = null; o.label = null; }
  else {
    if (!o.ref) { o.ref = Object.assign({ id: HM.uid(), label, source: o.point ? "point" : source, runs: o.runs }, o.point || {}, exportStats(o.st)); f.objects.push(o.ref); }
    else if (source === "manual" && o.ref.source === "similar") o.ref.source = "manual";
    o.ref.label = label; o.label = label;
  }
  F.suggest.delete(i); paintObj(i); HM.touch();
}
function addPoint(p, label) {
  const st = fs().settings, { w, h: hh } = F.px, rad = Math.max(3, Math.round(st.context_px / 2)), runs = [];
  const cx = Math.round(p[0]), cy = Math.round(p[1]);
  for (let y = Math.max(0, cy - rad); y <= Math.min(hh - 1, cy + rad); y++) { const dx = Math.floor(Math.sqrt(rad * rad - (y - cy) ** 2));
    let run = null;
    for (let x = Math.max(0, cx - dx); x <= Math.min(w - 1, cx + dx); x++) {
      if (F.lab[y * w + x] < 0) { if (run) run[2] = x; else run = [y, x, x]; } else if (run) { runs.push(...run); run = null; } }
    if (run) runs.push(...run); }
  if (!runs.length) { HM.status("No free pixels"); return; }
  const i = F.objs.length; paintRuns(F.lab, runs, i);
  const o = { ref: null, label: null, runs, point: { point: [cx, cy], radius_px: rad } }; o.st = stats(runs); F.objs.push(o);
  standardize(F.objs); setLabel(i, label, "point");
}
// For each unlabelled candidate, the nearest labelled example (up to 300 per class) in standardised feature space.
function findSimilar() {
  if (!F.lab) return;
  const ex = new Map();
  F.objs.forEach((o) => { if (o.label) { if (!ex.has(o.label)) ex.set(o.label, []); ex.get(o.label).push(o.v); } });
  for (const [k, list] of ex) if (list.length > 300) ex.set(k, list.slice(-300));
  if (![...ex.keys()].some((k) => k !== NONE)) { HM.status("No labelled examples"); return; }
  const thr = fs().settings.match, D = F.objs[0].v.length; F.suggest = new Map();
  F.objs.forEach((o, i) => { if (o.label) return; let best = Infinity, cls = null;
    for (const [k, list] of ex) for (const v of list) { let s = 0; for (let j = 0; j < D; j++) s += (o.v[j] - v[j]) ** 2;
      const d = Math.sqrt(s / D); if (d < best) { best = d; cls = k; } }
    if (cls && cls !== NONE && best < thr) F.suggest.set(i, cls); });
  paintAll(); renderPanel(); HM.status(`${F.suggest.size} suggested`);
}
function acceptSuggestions() { const n = F.suggest.size; for (const [i, c] of [...F.suggest]) setLabel(i, c, "similar");
  F.suggest.clear(); paintAll(); renderPanel(); HM.status(`${n} accepted`); }

// ---------- cartoon ----------
function colourOf(i) { const o = F.objs[i]; if (o.label) return C_LAB[o.label]; const s = F.suggest.get(i); return s ? C_SUG[s] : C_CAND; }
function paintAll() {
  if (!F.lab) return; const { w, h: hh, N } = F.px;
  if (!F.rimg) { F.raster = document.createElement("canvas"); F.raster.width = w; F.raster.height = hh;
    F.rctx = F.raster.getContext("2d"); F.rimg = F.rctx.createImageData(w, hh); }
  const d = F.rimg.data, lab = F.lab, R = F.region, T = F.tissue, cols = F.objs.map((_, i) => colourOf(i));
  for (let i = 0, j = 0; i < N; i++, j += 4) { const k = lab[i], c = k >= 0 ? cols[k] : R[i] ? C_TIS : T[i] ? C_OFF : C_OUT;
    d[j] = c[0]; d[j + 1] = c[1]; d[j + 2] = c[2]; d[j + 3] = 255; }
  F.rctx.putImageData(F.rimg, 0, 0); drawCartoon();
}
function paintObj(i) {
  if (!F.rimg) return; const o = F.objs[i], c = colourOf(i), d = F.rimg.data, w = F.px.w, r = o.runs;
  for (let k = 0; k < r.length; k += 3) for (let x = r[k + 1]; x <= r[k + 2]; x++) { const j = 4 * (r[k] * w + x); d[j] = c[0]; d[j + 1] = c[1]; d[j + 2] = c[2]; }
  const b = o.st.bbox; F.rctx.putImageData(F.rimg, 0, 0, b[0], b[1], b[2] - b[0] + 1, b[3] - b[1] + 1); drawCartoon();
}
// Same mapping as the SVG viewBox (xMidYMid meet), so both panes show the same part of the image.
function geom(cv) { const r = cv.getBoundingClientRect(), v = S.view, sc = Math.min(r.width / v.w, r.height / v.h);
  return { r, sc, ox: (r.width - v.w * sc) / 2 - v.x * sc, oy: (r.height - v.h * sc) / 2 - v.y * sc }; }
const toImgC = (ev) => { const g = geom($("cartoon")); return [(ev.clientX - g.r.left - g.ox) / g.sc, (ev.clientY - g.r.top - g.oy) / g.sc]; };
const cssColor = (v) => { const m = /--[\w-]+/.exec(v); return m ? getComputedStyle(document.documentElement).getPropertyValue(m[0]).trim() : v; };
function drawCartoon() {
  if (!F.on) return; const cv = $("cartoon"), g = geom(cv), dpr = window.devicePixelRatio || 1;
  const W = Math.round(g.r.width * dpr), H = Math.round(g.r.height * dpr); if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
  const x = cv.getContext("2d"); x.setTransform(1, 0, 0, 1, 0, 0); x.fillStyle = "#2a2a2e"; x.fillRect(0, 0, W, H);
  if (!S.doc) return;
  if (!(F.raster && current())) { x.fillStyle = "#bbb"; x.font = `${13 * dpr}px system-ui, sans-serif`;
    x.fillText(F.busy ? "Finding features…" : "Not detected", 14 * dpr, 24 * dpr); return; }
  x.setTransform(g.sc * dpr, 0, 0, g.sc * dpr, g.ox * dpr, g.oy * dpr);
  x.imageSmoothingEnabled = g.sc < 1; x.drawImage(F.raster, 0, 0);
  x.lineWidth = 1.6 / g.sc; x.lineJoin = "round";
  for (const c of S.doc.curves) { x.strokeStyle = cssColor(HM.KINDS[c.kind].color);
    for (const sg of c.segments) { if (sg.length < 2) continue; x.beginPath(); sg.forEach((p, i) => (i ? x.lineTo(p[0], p[1]) : x.moveTo(p[0], p[1])));
      if (HM.KINDS[c.kind].closed) x.closePath(); x.stroke(); } }
  if (F.hover >= 0 && F.objs[F.hover]) { const b = F.objs[F.hover].st.bbox, pad = 3 / g.sc;
    x.strokeStyle = "#facc15"; x.lineWidth = 2 / g.sc; x.strokeRect(b[0] - pad, b[1] - pad, b[2] - b[0] + 1 + 2 * pad, b[3] - b[1] + 1 + 2 * pad); }
}
function objAt(p) { if (!F.lab || !current()) return -1; const x = Math.floor(p[0]), y = Math.floor(p[1]);
  return x < 0 || y < 0 || x >= F.px.w || y >= F.px.h ? -1 : F.lab[y * F.px.w + x]; }
function setHover(i) { if (i === F.hover) return; F.hover = i; HM.renderAll(); }
function toggleCartoon(on = !F.on) {
  F.on = on; $("cartoon").hidden = !on; $("btnCartoon").classList.toggle("active", on);
  window.dispatchEvent(new Event("resize"));
  if (on && S.doc && !F.lab) detect(); else drawCartoon();
}

// ---------- popover for one feature ----------
const fmtArea = (px) => { const k = HM.pxPerUnit(); if (!k) return `${px} px²`;
  const u = S.doc.calibration.unit, a = px / k / k, um2 = u === "um" ? a : u === "mm" ? a * 1e6 : a * 1e8;
  return um2 < 1e5 ? `${Math.round(um2).toLocaleString()} µm²` : `${(um2 / 1e6).toFixed(3)} mm²`; };
function openFeature(i, p, ev) {
  const anchor = { getBoundingClientRect: () => ({ left: ev.clientX, right: ev.clientX, top: ev.clientY, bottom: ev.clientY }) };
  HM.showPop(anchor, () => renderFeature(i, p), "feat:" + i + ":" + Math.round(p[0]) + "," + Math.round(p[1]));
}
function renderFeature(i, p) {
  const m = $("pop"); m.innerHTML = "";
  const btn = (k, lab, n, on, f) => h("button", { class: on ? "on" : "", onclick: f },
    h("span", { class: "sw", style: `background:${HEX[k]}` }), lab, n ? h("kbd", {}, n) : null);
  if (i < 0) {
    m.append(h("div", { class: "tm-title" }, "No candidate here"),
      h("div", { class: "muted" }, "Add point:"));
    const s = h("div", { class: "seg col" });
    CLASSES.forEach(([k, lab]) => s.append(btn(k, lab, null, false, () => { addPoint(p, k); renderPanel(); HM.closePop(); })));
    m.append(s, h("div", { class: "tm-foot" }, h("button", { class: "primary", onclick: HM.closePop }, "Done"))); return;
  }
  const o = F.objs[i], st = o.st, sug = F.suggest.get(i);
  const shape = st.elong > 0.6 ? (st.radial == null ? "elongate" : st.radial > 0.7 ? "elongate, radial" : st.radial < 0.3 ? "elongate, circumferential" : "elongate, oblique") : "roundish";
  m.append(h("div", { class: "tm-title" }, o.label ? NAME[o.label] : sug ? `Suggested: ${NAME[sug]}` : "Unlabelled candidate"),
    h("div", { class: "muted" }, `${fmtArea(st.n)} · ${shape}${o.ref ? " · " + o.ref.source : ""}`));
  const s = h("div", { class: "seg col" }); const re = () => { renderFeature(i, p); renderPanel(); };
  CLASSES.forEach(([k, lab], n) => s.append(btn(k, lab, n + 1, o.label === k, () => { setLabel(i, k); re(); })));
  s.append(btn(NONE, NAME[NONE], "X", o.label === NONE, () => { setLabel(i, NONE); re(); }));
  m.append(s, h("div", { class: "tm-foot" },
    sug ? h("button", { onclick: () => { setLabel(i, sug, "similar"); re(); } }, "Accept suggestion") : null,
    o.label ? h("button", { onclick: () => { setLabel(i, null); re(); } }, "Clear") : null,
    h("button", { onclick: () => { HM.closePop(); findSimilar(); } }, "Find similar"),
    h("button", { class: "primary", onclick: HM.closePop }, "Done")));
}

// ---------- panel ----------
function renderPanel() {
  const box = $("featPanel"); if (!box) return; box.innerHTML = ""; if (!S.doc) return;
  const f = fs(), st = f.settings, re = () => { HM.touch(); schedule(); renderPanel(); };
  const num = (k, title, min, max) => h("input", { type: "number", class: "num", min, max, step: 1, value: st[k], title,
    onchange: (e) => { const v = parseFloat(e.target.value); if (isFinite(v)) { st[k] = Math.min(max, Math.max(min, v)); re(); } } });
  const sens = h("input", { type: "range", min: 0, max: 100, value: Math.round((THR[0] - st.threshold) / (THR[0] - THR[1]) * 100),
    oninput: (e) => { st.threshold = +(THR[0] - (THR[0] - THR[1]) * e.target.value / 100).toFixed(2); HM.touch(); schedule(); } });
  box.append(h("div", { class: "grid2" },
    h("label", {}, "Light"), HM.seg(["ppl", "xpl"], st.light, (v) => { st.light = v || st.light; HM.touch(); renderPanel(); }, (v) => (v === "ppl" ? "plain" : "XPL")),
    h("label", {}, "Look for"),
    HM.seg(["dark", "bright", "either"], st.polarity, (v) => { st.polarity = v || st.polarity; re(); }, (v) => ({ dark: "darker", bright: "brighter", either: "either" })[v]),
    h("label", {}, "Sensitivity"), sens,
    h("label", {}, "Window px"),
    h("span", {}, num("window_px", "", 3, 80), h("span", { class: "muted" }, " min area "), num("min_area_px", "", 1, 5000))));
  const where = { areas: "areas", cortex: "cortex", tissue: "all tissue" }[F.regionKind];
  const objs = current() ? F.objs : [], reg = F.regionPx || 1;
  const cand = objs.filter((o) => !o.label), candA = cand.reduce((a, o) => a + o.st.n, 0);
  box.append(h("div", { class: "box" }, F.busy ? "Finding features…" : !F.lab || !current() ? "—" :
    `${cand.length} unlabelled · ${(100 * candA / reg).toFixed(1)}% · ${where}`));
  const ul = h("ul", { class: "legend" });
  CLASSES.forEach(([k, lab], n) => { const mine = f.objects.filter((o) => o.label === k), a = mine.reduce((s, o) => s + o.area_px, 0);
    const sug = [...F.suggest.values()].filter((c) => c === k).length;
    ul.append(h("li", {}, h("kbd", {}, String(n + 1)), h("span", { class: "sw", style: `background:${HEX[k]}` }), h("span", { class: "kind" }, lab),
      h("span", { class: "detail" }, (mine.length ? `${mine.length}` + (F.regionPx ? ` · ${(100 * a / reg).toFixed(1)}%` : "") : "") + (sug ? ` +${sug}?` : "")))); });
  const nn = f.objects.filter((o) => o.label === NONE).length;
  ul.append(h("li", {}, h("kbd", {}, "X"), h("span", { class: "sw", style: `background:${HEX[NONE]}` }), h("span", { class: "kind" }, NAME[NONE]),
    h("span", { class: "detail" }, nn ? String(nn) : "")));
  box.append(ul);
  box.append(h("div", { class: "tools" },
    h("button", { onclick: findSimilar, title: "F" }, "Find similar"),
    h("input", { type: "range", min: 30, max: 300, value: Math.round(st.match * 100), title: "Match",
      oninput: (e) => { st.match = e.target.value / 100; }, onchange: () => { if (F.suggest.size) findSimilar(); } }),
    h("button", { class: F.suggest.size ? "primary" : "", disabled: !F.suggest.size, onclick: acceptSuggestions, title: "A" },
      `Accept ${F.suggest.size || ""}`.trim())));
}

// ---------- mask (export and reload) ----------
// Labelled objects get values 1..n in doc.features.objects order; unlabelled candidates follow.
function beforeExport(d, base) {
  const f = fs(d); f.objects.forEach((o, i) => { o.mask_value = i + 1; });
  const cands = current() ? F.objs.filter((o) => !o.ref).length : 0;
  f.mask = f.objects.length || cands ? { file: base + ".mask.png", encoding: MASK_ENCODING, n_labelled: f.objects.length, n_unlabelled: cands } : null;
}
async function maskFile(d, base) {
  const f = fs(d); if (!f.mask) return null;
  const w = S.img.w, hh = S.img.h, cv = document.createElement("canvas"); cv.width = w; cv.height = hh;
  const g = cv.getContext("2d"), im = g.createImageData(w, hh), a = im.data;
  for (let i = 3; i < a.length; i += 4) a[i] = 255;
  const put = (r, v) => { for (let k = 0; k < r.length; k += 3) for (let x = r[k + 1]; x <= r[k + 2]; x++) { const j = 4 * (r[k] * w + x);
    a[j] = v & 255; a[j + 1] = (v >> 8) & 255; a[j + 2] = (v >> 16) & 255; } };
  let v = f.objects.length;
  if (current()) for (const o of F.objs) if (!o.ref) put(o.runs, ++v);
  f.objects.forEach((o) => put(o.runs, o.mask_value));
  g.putImageData(im, 0, 0);
  const blob = await new Promise((r) => cv.toBlob(r, "image/png"));
  return { name: f.mask.file, blob };
}
async function loadMask(file) {
  if (!(await ensurePixels())) return; prepare(fs().settings);
  const url = URL.createObjectURL(file), im = new Image(); im.src = url;
  try { await im.decode(); } finally { URL.revokeObjectURL(url); }
  const { w, h: hh, N } = F.px;
  if (im.naturalWidth !== w || im.naturalHeight !== hh) { HM.status("Mask size differs; ignored", true); return detect(); }
  const cv = document.createElement("canvas"); cv.width = w; cv.height = hh;
  const g = cv.getContext("2d", { willReadFrequently: true }); g.drawImage(im, 0, 0); const a = g.getImageData(0, 0, w, hh).data;
  const byV = new Map();
  for (let y = 0; y < hh; y++) for (let x = 0; x < w;) { const j = 4 * (y * w + x), v = a[j] + 256 * a[j + 1] + 65536 * a[j + 2];
    if (!v) { x++; continue; } const x0 = x; while (x < w && (() => { const k = 4 * (y * w + x); return a[k] + 256 * a[k + 1] + 65536 * a[k + 2] === v; })()) x++;
    if (!byV.has(v)) byV.set(v, []); byV.get(v).push(y, x0, x - 1); }
  const f = fs(), refs = new Map(f.objects.map((o) => [o.mask_value, o])), lab = new Int32Array(N).fill(-1), objs = [];
  for (const [v, runs] of byV) { const ref = refs.get(v) || null; if (ref) ref.runs = runs;
    paintRuns(lab, runs, objs.length); objs.push({ ref, label: ref ? ref.label : null, runs }); }
  for (const o of f.objects) if (!byV.has(o.mask_value)) { paintRuns(lab, o.runs, objs.length); objs.push({ ref: o, label: o.label, runs: o.runs }); }
  F.fromMask = true; finish(lab, objs);
  HM.status(`Mask: ${f.objects.length} labelled, ${objs.length - f.objects.length} unlabelled`);
}

// ---------- derived numbers for the export ----------
// Zone of a point = the nearest boundary outward along the ray from the centre (descriptors live on that boundary).
function zoneFinder(d) {
  const z = d.centre; if (!z) return () => null;
  const bs = d.curves.filter((c) => HM.KINDS[c.kind].boundary && !c.branch_of && c.segments.some((s) => s.length > 1));
  const prof = bs.map((c) => { const acc = Array.from({ length: 72 }, () => [0, 0]);
    for (const q of c.segments.flatMap((s) => HM.resample(s, 2))) { const b = Math.floor(HM.angle(q, z) / 5) % 72; acc[b][0] += HM.dist(q, z); acc[b][1]++; }
    return acc.map(([s, n]) => (n ? s / n : null)); });
  return (x, y) => { const q = [x, y], b = Math.floor(HM.angle(q, z) / 5) % 72, r = HM.dist(q, z); let best = null, br = Infinity;
    bs.forEach((c, i) => { const pr = prof[i][b]; if (pr != null && pr >= r && pr < br) { br = pr; best = c.id; } }); return best; };
}
function deriveFeatures(d) {
  const f = d.features; if (!f || (!f.objects.length && !(d === S.doc && current() && F.lab))) return null;
  const k = HM.pxPerUnit(d), unit = k ? d.calibration.unit : "px", ar = (px) => (k ? px / k / k : px);
  const num = (v) => (v == null || !isFinite(v) ? null : +v.toPrecision(6)), zoneOf = zoneFinder(d), zs = new Map();
  const Z = (id) => { if (!zs.has(id)) zs.set(id, { px: 0, by: {} }); return zs.get(id); };
  const haveRegion = d === S.doc && current() && F.region;
  if (haveRegion) { const { w, h: hh } = F.px, s = 3;
    for (let y = 1; y < hh; y += s) for (let x = 1; x < w; x += s) if (F.region[y * w + x]) Z(zoneOf(x, y)).px += s * s; }
  for (const o of f.objects) { if (o.label === NONE) continue; const b = Z(zoneOf(o.centroid[0], o.centroid[1])).by;
    const e = b[o.label] || (b[o.label] = { n: 0, area_px: 0 }); e.n++; e.area_px += o.area_px; }
  const cls = (by, px) => Object.fromEntries(Object.entries(by).map(([l, e]) => [l, { n: e.n, area: num(ar(e.area_px)), area_fraction: px ? num(e.area_px / px) : null }]));
  const all = {}; for (const z of zs.values()) for (const [l, e] of Object.entries(z.by)) { const t = all[l] || (all[l] = { n: 0, area_px: 0 }); t.n += e.n; t.area_px += e.area_px; }
  const name = (id) => { const c = id && d.curves.find((x) => x.id === id); return c ? HM.curveName(c) : null; };
  return { note: "Region areas are sampled on a 3-px grid. Zone = nearest boundary outward from the centre.",
    unit2: unit + "2", region: haveRegion ? F.regionKind : null, region_area: haveRegion ? num(ar(F.regionPx)) : null,
    by_class: cls(all, haveRegion ? F.regionPx : 0),
    by_zone: [...zs.entries()].map(([id, z]) => ({ zone_outer_id: id, zone_outer: name(id), region_area: haveRegion ? num(ar(z.px)) : null, classes: cls(z.by, z.px) })) };
}

HM.features = { F, detect, setLabel, findSimilar, acceptSuggestions };   // for scripted tests

// ---------- wiring ----------
HM.hooks.open.push(() => { reset(); renderPanel(); if (S.doc && (F.on || fs().objects.length)) detect(); else drawCartoon(); });
HM.hooks.loaded.push((mask) => { reset(); renderPanel(); if (mask) loadMask(mask); else if (F.on || fs().objects.length) detect(); });
HM.hooks.render.push(() => {
  if (F.hover >= 0 && F.objs[F.hover] && current()) { const b = F.objs[F.hover].st.bbox, k = HM.upp(), pad = 3 * k;
    HM.el("rect", { x: b[0] - pad, y: b[1] - pad, width: b[2] - b[0] + 1 + 2 * pad, height: b[3] - b[1] + 1 + 2 * pad, fill: "none",
      stroke: "var(--sel)", "stroke-width": 2 * k, "pointer-events": "none" }, HM.layerVerts); }
  drawCartoon();
});
HM.hooks.label.push((p, ev) => { if (!F.lab) { HM.status("Not detected yet"); return; } const i = objAt(p); openFeature(i, p, ev); });
HM.hooks.beforeExport.push(beforeExport);
HM.hooks.exportFiles.push(maskFile);
HM.hooks.derive.push((d, out) => { const x = deriveFeatures(d); if (x) out.features = x; });

const cv = $("cartoon"); let down = null;
cv.addEventListener("contextmenu", (e) => e.preventDefault());
cv.addEventListener("pointerdown", (e) => { if (!S.doc) return; HM.touch(); down = { x: e.clientX, y: e.clientY, v: Object.assign({}, S.view), moved: false }; cv.setPointerCapture(e.pointerId); });
cv.addEventListener("pointermove", (e) => {
  if (down) { const dx = e.clientX - down.x, dy = e.clientY - down.y; if (Math.hypot(dx, dy) > 4) down.moved = true;
    if (down.moved) { const sc = geom(cv).sc; S.view.x = down.v.x - dx / sc; S.view.y = down.v.y - dy / sc; HM.setView(); } return; }
  setHover(objAt(toImgC(e)));
});
cv.addEventListener("pointerup", (e) => { const d = down; down = null; if (!d || d.moved || !S.doc) return;
  if (!F.lab) { detect(); return; } const p = toImgC(e); openFeature(objAt(p), p, e); });
cv.addEventListener("pointerleave", () => setHover(-1));
cv.addEventListener("wheel", (e) => { if (!S.doc) return; e.preventDefault(); HM.zoomAt(toImgC(e), Math.exp(e.deltaY * (e.ctrlKey ? 0.01 : 0.0015))); }, { passive: false });
HM.svg.addEventListener("pointermove", (e) => { if (F.lab && (S.tool === "label" || F.on) && !S.drag) setHover(objAt(HM.toImg(e))); });
HM.svg.addEventListener("pointerleave", () => setHover(-1));

$("btnDetect").addEventListener("click", detect);
$("btnCartoon").addEventListener("click", () => toggleCartoon());
document.addEventListener("keydown", (ev) => {
  if (/^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement.tagName) || ev.metaKey || ev.ctrlKey || ev.altKey || !S.doc) return;
  const k = ev.key.toLowerCase();
  if (k === "v") return toggleCartoon();
  if (k === "l") return HM.setTool(S.tool === "label" ? "pan" : "label");
  if (k === "f") return findSimilar();
  if (k === "a" && F.suggest.size) return acceptSuggestions();
  if (F.hover < 0 || !F.lab) return;
  const n = parseInt(k, 10);
  if (n >= 1 && n <= CLASSES.length) setLabel(F.hover, CLASSES[n - 1][0]);
  else if (k === "x") setLabel(F.hover, NONE);
  else if (k === "u") setLabel(F.hover, null);
  else return;
  renderPanel();
});
})();
