/* Histology Measure: guided observations.
 * Students answer plain questions about what they see. A versioned rule table turns the answers into the histological
 * descriptors the importer and the model use (matrix, vasc_density, vasc_orient, remodelling, tissue). The export keeps
 * the raw answers (`obs`), the derived terms and the rule version, so the rules can be revised and re-applied at import
 * without measuring again. Loaded before app.js; no dependencies.
 */
(function () {
"use strict";
const RULES_VERSION = "obs-1";

// ---------- icons: a small piece of cortex, outer surface at the top ----------
const box = (inner) => `<svg viewBox="0 0 40 32" width="40" height="32" aria-hidden="true"><rect x="1" y="4" width="38" height="27" rx="3" class="ic-bone"/>` +
  `<path d="M1 4.5 H39" class="ic-surf"/>${inner}</svg>`;
const dots = (pts, r = 1.6) => pts.map(([x, y]) => `<circle cx="${x}" cy="${y}" r="${r}" class="ic-canal"/>`).join("");
const lines = (ls) => ls.map((d) => `<path d="${d}" class="ic-line"/>`).join("");
const ICONS = {
  amt_none: box(""),
  amt_few: box(dots([[9, 12], [27, 20], [16, 26]])),
  amt_many: box(dots([[7, 10], [17, 11], [28, 10], [12, 18], [23, 18], [33, 19], [7, 26], [18, 26], [29, 26]])),
  amt_packed: box(dots([[6, 9], [13, 10], [20, 9], [27, 10], [34, 9], [9, 15], [16, 16], [23, 15], [30, 16], [6, 21], [13, 22], [20, 21],
    [27, 22], [34, 21], [9, 27], [16, 27], [23, 27], [30, 27]], 2.2)),
  pat_dots: box(dots([[8, 11], [18, 12], [29, 10], [12, 19], [24, 20], [33, 18], [8, 26], [19, 27], [30, 26]])),
  pat_circ: box(lines(["M5 11 H17 M21 11 H35", "M4 18 H14 M18 18 H36", "M6 25 H22 M26 25 H35"])),
  pat_radial: box(lines(["M8 8 V16 M8 20 V28", "M16 9 V27", "M24 8 V14 M24 18 V28", "M32 9 V27"])),
  pat_grid: box(lines(["M4 11 H36", "M4 18 H36", "M4 25 H36", "M10 11 V18", "M26 11 V18", "M17 18 V25", "M32 18 V25"])),
  pat_net: box(lines(["M5 10 L12 15 L9 22 L16 27", "M12 15 L21 12 L27 18 L22 25", "M27 18 L35 13", "M21 12 L24 7", "M9 22 L22 25 L33 27"])),
  pat_mixed: box(dots([[8, 11], [14, 25], [31, 25]]) + lines(["M18 11 H35", "M22 16 V24", "M4 18 H14"])),
  ost_yes: box(`<circle cx="12" cy="17" r="5" class="ic-halo"/><circle cx="12" cy="17" r="1.8" class="ic-canal"/>` +
    `<circle cx="28" cy="20" r="5" class="ic-halo"/><circle cx="28" cy="20" r="1.8" class="ic-canal"/>`),
  ost_no: box(dots([[12, 17], [28, 20]], 1.8)),
  sec_none: box(dots([[10, 14], [27, 22]])),
  sec_few: box(`<circle cx="15" cy="17" r="6" class="ic-sec"/><circle cx="15" cy="17" r="1.8" class="ic-canal"/>` + dots([[30, 12], [31, 25]])),
  sec_many: box([[9, 13], [21, 21], [32, 12], [11, 25]].map(([x, y]) => `<circle cx="${x}" cy="${y}" r="5" class="ic-sec"/><circle cx="${x}" cy="${y}" r="1.5" class="ic-canal"/>`).join("")),
  sec_most: box([[7, 10], [17, 10], [28, 11], [36, 11], [11, 19], [22, 20], [32, 21], [6, 27], [17, 28], [28, 28]].map(([x, y]) =>
    `<circle cx="${x}" cy="${y}" r="4.6" class="ic-sec"/><circle cx="${x}" cy="${y}" r="1.3" class="ic-canal"/>`).join("")),
  mx_woven: box([...Array(42)].map((_, i) => `<circle cx="${3 + (i * 7.3) % 34}" cy="${7 + (i * 4.1) % 23}" r="0.9" class="ic-speck"/>`).join("")),
  mx_pf: box(`<path d="M3 9 H37 M3 13 H37 M3 17 H37 M3 21 H37 M3 25 H37 M3 29 H37" class="ic-fibre"/>`),
  mx_lam: box(`<rect x="1" y="8" width="38" height="3" class="ic-bright"/><rect x="1" y="14" width="38" height="3" class="ic-bright"/>` +
    `<rect x="1" y="20" width="38" height="3" class="ic-bright"/><rect x="1" y="26" width="38" height="3" class="ic-bright"/>`),
  mx_fl: box([...Array(30)].map((_, i) => `<circle cx="${3 + (i * 7.3) % 34}" cy="${7 + (i * 4.1) % 23}" r="0.9" class="ic-speck"/>`).join("") +
    `<circle cx="12" cy="17" r="5" class="ic-ring"/><circle cx="12" cy="17" r="1.8" class="ic-canal"/><circle cx="29" cy="21" r="5" class="ic-ring"/><circle cx="29" cy="21" r="1.8" class="ic-canal"/>`),
  lac_plump: box([[8, 11, 20], [18, 14, 70], [29, 10, 130], [12, 22, 160], [25, 24, 40], [33, 19, 100]].map(([x, y, a]) =>
    `<ellipse cx="${x}" cy="${y}" rx="2.4" ry="1.7" transform="rotate(${a} ${x} ${y})" class="ic-canal"/>`).join("")),
  lac_flat: box([[8, 11], [19, 11], [30, 11], [13, 18], [24, 18], [34, 18], [8, 25], [19, 25], [30, 25]].map(([x, y]) =>
    `<ellipse cx="${x}" cy="${y}" rx="3.2" ry="0.9" class="ic-canal"/>`).join("")),
  efs_yes: box(lines(["M2 7 H38", "M2 9.5 H38", "M2 12 H38"]) + dots([[10, 20], [25, 24], [33, 18]])),
  efs_no: box(dots([[8, 9], [20, 8], [31, 10], [10, 20], [25, 24], [33, 18]])),
};

// ---------- questions ----------
// opts: [value, label, icon]. `when(obs, ctx)` hides a question that does not apply.
const ZONE_Q = [
  { key: "canal_amount", q: "How much of this bone is canals (holes and channels)?",
    opts: [["none", "None, or almost none", "amt_none"], ["few", "A few, scattered", "amt_few"],
           ["many", "Many, spread through it", "amt_many"], ["packed", "Packed: almost spongy", "amt_packed"]] },
  { key: "canal_pattern", q: "What do most canals look like?",
    when: (o) => o.canal_amount && o.canal_amount !== "none",
    opts: [["dots", "Round dots or short ovals", "pat_dots"], ["circ", "Long lines running around the bone", "pat_circ"],
           ["radial", "Lines running outward, like spokes", "pat_radial"], ["grid", "Lines around, joined by spokes", "pat_grid"],
           ["net", "An irregular net", "pat_net"], ["mixed", "A mix", "pat_mixed"]] },
  { key: "prim_osteons", q: "Do the canals have a halo of layered bone around them, with no sharp outer edge?",
    when: (o) => o.canal_amount && o.canal_amount !== "none",
    opts: [["yes", "Yes, most do", "ost_yes"], ["no", "No, plain holes", "ost_no"]] },
  { key: "sec_osteons", q: "Are there round osteons with a sharp outer ring that cut across the tissue around them?",
    opts: [["none", "None", "sec_none"], ["few", "A few", "sec_few"], ["many", "Many", "sec_many"], ["most", "The bone is mostly made of them", "sec_most"]] },
  { key: "matrix_xpl", q: "Under crossed polars, the bone between the canals is…", when: (o, c) => c.polarized,
    opts: [["woven", "Dark and grainy (speckled)", "mx_woven"], ["pf", "Bright, fibres running one way", "mx_pf"],
           ["lam", "Fine bright and dark layers", "mx_lam"], ["fl", "Grainy, with bright rings around canals", "mx_fl"]] },
  { key: "lacuna_shape", q: "The tiny dark specks (osteocyte lacunae) are…", when: (o, c) => !c.polarized,
    opts: [["plump", "Plump, round, pointing every way", "lac_plump"], ["flat", "Flat and lined up in rows", "lac_flat"],
           ["mixed", "Both kinds", null]] },
];

const BONE_Q = [
  { key: "efs", q: "Just under the outer surface, are several growth lines packed tightly together, with few or no canals between them?",
    opts: [["present", "Yes, packed lines at the surface", "efs_yes"], ["absent", "No: ordinary tissue reaches the surface", "efs_no"]],
    reasons: ["outer surface broken or eroded", "outer surface not in the image", "image too coarse", "other"] },
  { key: "sec_extent", q: "How far out from the marrow cavity do secondary osteons reach?",
    opts: [["none", "No secondary osteons"], ["inner", "Inner cortex only"], ["mid", "Into the middle cortex"], ["outer", "Out to the surface"]],
    reasons: ["inner cortex not in the image", "image too coarse", "other"] },
  { key: "icl", q: "Is the marrow cavity lined by a smooth band of layered bone, with a sharp line on its outer side?",
    opts: [["present", "Yes"], ["absent", "No"]], reasons: ["marrow cavity not in the image", "image too coarse", "other"] },
  { key: "medullary_bone", q: "Inside the marrow cavity, is there a separate lining of spongy, woven-looking bone?",
    opts: [["present", "Yes"], ["absent", "No"]], reasons: ["marrow cavity not in the image", "image too coarse", "other"] },
  { key: "kastschenko_line", q: "Is there a wavy cement line where cortex meets the inner (endosteal) bone?", more: true,
    opts: [["present", "Yes"], ["absent", "No"]], reasons: ["not in the image", "other"] },
  { key: "image_quality", q: "How readable is this image overall?",
    opts: [["good", "Good"], ["fair", "Fair"], ["poor", "Poor"]] },
];

// ---------- interpretation ----------
const DENS = { none: "none", few: "low", many: "moderate", packed: "high" };
const ORIENT = { dots: "longitudinal", circ: "circumferential", radial: "radial", grid: "plexiform", net: "reticular", mixed: "mixed" };
const REMOD = { none: "none", few: "scattered", many: "dense", most: "dense" };
const MX_XPL = { woven: "woven", pf: "parallel_fibred", lam: "lamellar", fl: "woven" };
const MX_LAC = { plump: "woven", flat: "organised", mixed: "mixed" };
const UNCLEAR = "unclear";

// obs -> {matrix, vasc_density, vasc_orient, remodelling, prim_osteons, tissue, basis_light}; null where unanswered.
function interpretZone(o, ctx) {
  o = o || {};
  const out = { rules_version: RULES_VERSION };
  const a = o.canal_amount;
  out.vasc_density = a === UNCLEAR ? UNCLEAR : DENS[a] || null;
  if (a === "none") out.vasc_orient = null;
  else if (o.canal_pattern === UNCLEAR) out.vasc_orient = UNCLEAR;
  else if (o.canal_pattern === "circ") out.vasc_orient = a === "many" || a === "packed" ? "laminar" : "circumferential";
  else out.vasc_orient = ORIENT[o.canal_pattern] || null;
  out.remodelling = o.sec_osteons === UNCLEAR ? UNCLEAR : REMOD[o.sec_osteons] || null;
  let mx = null, basis = null;
  if (ctx.polarized && o.matrix_xpl) { mx = o.matrix_xpl === UNCLEAR ? UNCLEAR : MX_XPL[o.matrix_xpl]; basis = "crossed_polars"; }
  else if (!ctx.polarized && o.lacuna_shape) { mx = o.lacuna_shape === UNCLEAR ? UNCLEAR : MX_LAC[o.lacuna_shape]; basis = "plain_light_lacunae"; }
  out.matrix = mx; out.basis_light = basis;
  const po = o.matrix_xpl === "fl" || o.prim_osteons === "yes" ? true : o.prim_osteons === "no" ? false : null;
  out.prim_osteons = a === "none" ? false : po;
  out.tissue = tissue(out, o);
  return out;
}
function tissue(t, o) {
  if (o.sec_osteons === "most") return "secondary_haversian";
  const m = t.matrix, po = t.prim_osteons, v = t.vasc_density;
  if (!m || m === UNCLEAR) return null;
  if (m === "woven") return po ? "fibrolamellar" : po === false ? "woven_fibred" : null;
  if (m === "parallel_fibred") return po ? "parallel_fibred_with_primary_osteons" : "parallel_fibred";
  if (m === "lamellar") return v === "none" || v === "low" ? "lamellar" : "lamellar_vascular";
  if (m === "organised") return "parallel_fibred_or_lamellar";
  if (m === "mixed") return "mixed";
  return null;
}
const NICE = { secondary_haversian: "secondary (Haversian) bone", fibrolamellar: "fibrolamellar", woven_fibred: "woven-fibred",
  parallel_fibred_with_primary_osteons: "parallel-fibred with primary osteons", parallel_fibred: "parallel-fibred", lamellar: "lamellar",
  lamellar_vascular: "vascular lamellar", parallel_fibred_or_lamellar: "parallel-fibred or lamellar", mixed: "mixed matrix" };
// One line describing the derived terms (not shown to measurers; for scripts and checks).
function summary(t) {
  if (!t) return "";
  const parts = [];
  if (t.tissue) parts.push(NICE[t.tissue] || t.tissue.replace(/_/g, " "));
  else if (t.matrix && t.matrix !== UNCLEAR) parts.push(t.matrix.replace(/_/g, "-") + " matrix");
  if (t.vasc_orient && t.vasc_orient !== UNCLEAR) parts.push(t.vasc_orient + " canals");
  if (t.vasc_density) parts.push(t.vasc_density === "none" ? "avascular" : t.vasc_density === UNCLEAR ? "" : t.vasc_density + " vascularity");
  if (t.remodelling && t.remodelling !== "none" && t.remodelling !== UNCLEAR) parts.push(t.remodelling + " secondary osteons");
  return parts.filter(Boolean).join(" · ");
}
// Answered every question that applies?
function complete(o, ctx) {
  o = o || {};
  return ZONE_Q.every((q) => (q.when && !q.when(o, ctx)) || o[q.key] != null);
}

window.HG = { RULES_VERSION, ICONS, ZONE_Q, BONE_Q, interpretZone, summary, complete, UNCLEAR };
})();
