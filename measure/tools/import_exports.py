#!/usr/bin/env python3
"""Import Histology Measure exports (JSON, schema 0.1-0.4) into flat tables for analysis.

Everything is recomputed from the raw image-pixel points and the scale bar; the `derived` block in each export
is used only as a cross-check. All lengths are written in mm (um and cm exports are converted).

Usage:
    python3 measure/tools/import_exports.py EXPORTS [EXPORTS ...] [--out data_clean] [--prefix digitized_]
EXPORTS are export files or folders (searched recursively for *.json files that carry `curves`).

Writes <out>/<prefix>{sections,curves,zones,transects,lacunae,features,feature_zones}.csv:
  sections  one row per kept export: identity, image setup (view, light, other-light image), section metadata,
            flags, scores (with "can't tell" reasons), mark counts as an interval [count_lo, count_hi], outer-edge
            preservation (periosteum vs broken edge), skipped steps and their reasons, scale, import issues
  curves    one row per traced curve, innermost first (`radial_rank`), with length, radii, angular coverage and
            two circumference estimates (see README, "Import")
  zones     one row per zone between consecutive boundaries (or one "field" row for a close-up), with the
            student's raw answers (obs_*) and the descriptors derived from them by guide.js (rules_version)
  transects one row per transect crossing: order, kind, mark number, distance from the outer end, radius from
            the centre and width to the next crossing inward (mm)
  lacunae   one row per counting area: lacunae counted, area (mm2), density (per mm2)
  features  one row per labelled feature (canal, osteon, ...): label, centroid, area (mm2) and shape, recomputed from
            its pixel runs, and the zone it lies in (nearest boundary outward from the centre)
  feature_zones  per zone and label: count, area and area fraction, with the zone's searched area taken from the
            tool's `derived.features` block (it needs the image, so it cannot be recomputed here)

If the same measurer exported the same image more than once (same figure_id, measurer and repeat flag), only the
latest export is kept and the others are listed as superseded. Standard library only.
"""
import argparse, csv, json, math, os, sys
from collections import defaultdict

BOUNDARY = {"periosteum", "broken_edge", "medullary", "remodel_front", "mark", "label", "efs_inner"}
OUTER = {"periosteum", "broken_edge"}
OBS = ("canal_amount", "canal_pattern", "prim_osteons", "sec_osteons", "matrix_xpl", "lacuna_shape")
DESC = ("matrix", "vasc_density", "vasc_orient", "remodelling", "prim_osteons", "tissue", "basis_light", "rules_version")
TO_MM = {"mm": 1.0, "um": 1e-3, "cm": 10.0}
BIN = 5                      # degrees per angular bin, as in the tool
NB = 360 // BIN


# ---------- geometry (mirrors measure/app.js) ----------
def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def resample(seg, step=2.0):
    out = []
    for a, b in zip(seg, seg[1:]):
        n = max(1, math.ceil(dist(a, b) / step))
        out += [[a[0] + (b[0] - a[0]) * j / n, a[1] + (b[1] - a[1]) * j / n] for j in range(n)]
    if seg:
        out.append(seg[-1])
    return out


def sample(c, mx=60):
    p = [q for s in c["segments"] for q in resample(s)]
    if len(p) <= mx:
        return p
    k = len(p) / mx
    return [p[int(i * k)] for i in range(mx)]


def angle(q, z):
    return (math.degrees(math.atan2(q[1] - z[1], q[0] - z[0])) + 360) % 360


def dist_to_polyline(p, segs):
    best = math.inf
    for s in segs:
        for a, b in zip(s, s[1:]):
            dx, dy = b[0] - a[0], b[1] - a[1]
            L2 = dx * dx + dy * dy
            t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2)) if L2 else 0.0
            best = min(best, math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy))
    return best


def poly_area(pts):
    return abs(sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(pts, pts[1:] + pts[:1]))) / 2


def in_poly(p, pts):
    inside, j = False, len(pts) - 1
    for i in range(len(pts)):
        a, b = pts[i], pts[j]
        if (a[1] > p[1]) != (b[1] > p[1]) and p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]:
            inside = not inside
        j = i
    return inside


def binned_radius(c, z):
    """Mean radius (px) from the centre in each 5-degree bin; None where the curve was not traced."""
    acc = defaultdict(list)
    for s in c["segments"]:
        for q in resample(s):
            acc[int(angle(q, z) // BIN) % NB].append(dist(q, z))
    return [sum(acc[i]) / len(acc[i]) if acc[i] else None for i in range(NB)]


def arc(bins):
    """Coverage (deg) and occupied arc [start, end], same rule as the tool: the largest empty run is the gap."""
    n = sum(b is not None for b in bins)
    if n == NB:
        return 360, 0, 360
    gs, gl = 0, 0
    for i in range(NB):
        if bins[i] is None:
            L = 0
            while L < NB and bins[(i + L) % NB] is None:
                L += 1
            if L > gl:
                gs, gl = i, L
    return n * BIN, ((gs + gl) % NB) * BIN, (gs * BIN) % 360 or 360


def fill_circular(bins):
    """Linear interpolation of missing bins around the circle."""
    have = [i for i, b in enumerate(bins) if b is not None]
    if not have:
        return None
    out = list(bins)
    for i in range(NB):
        if out[i] is None:
            lo = max((h for h in have if h < i), default=have[-1] - NB)
            hi = min((h for h in have if h > i), default=have[0] + NB)
            t = (i - lo) / (hi - lo)
            out[i] = bins[lo % NB] * (1 - t) + bins[hi % NB] * t
    return out


def star_perimeter(r):
    pts = [(r[i] * math.cos(math.radians((i + .5) * BIN)), r[i] * math.sin(math.radians((i + .5) * BIN))) for i in range(NB)]
    return sum(dist(p, q) for p, q in zip(pts, pts[1:] + pts[:1]))


# ---------- reading ----------
def migrate(d):
    if d.get("schema_version") in ("0.2", "0.3", "0.4"):
        d.setdefault("features", {"settings": {}, "objects": []})
        return d, None
    if d.get("schema_version") != "0.1":
        return None, f"unknown schema_version {d.get('schema_version')!r}"
    for c in d.get("curves", []):
        if c.get("kind") == "mark":
            if c.get("mark_type") == "double_multiple":
                c["mark_type"], c["branch_of"] = "LAG", c.get("linked_to")
            c.pop("linked_to", None)
            c.setdefault("branch_of", None)
    d.setdefault("section", {})
    if (d.get("flags") or {}).pop("section_oblique", False):
        d["section"]["plane"] = "oblique"
    d.setdefault("lacunae", [])
    d["migrated_from"] = "0.1"
    d["features"] = {"settings": {}, "objects": []}
    return d, None


def find_exports(paths):
    for p in paths:
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                for f in sorted(files):
                    if f.endswith(".json") and f != "manifest.json" and f != "truth.json":
                        yield os.path.join(root, f)
        else:
            yield p


def scale(d):
    c = d.get("calibration") or {}
    if not (c.get("p1") and c.get("p2") and (c.get("length") or 0) > 0):
        return None, None
    u = c.get("unit", "mm")
    if u not in TO_MM:
        return None, u
    return dist(c["p1"], c["p2"]) / (c["length"] * TO_MM[u]), u      # px per mm


def r4(v):
    return None if v is None or not math.isfinite(v) else round(v, 4)


def segs_ok(c):
    return any(len(s) for s in c.get("segments", []))


def outer_segments(d):
    return [s for c in d["curves"] if c["kind"] in OUTER for s in c["segments"] if len(s) > 1]


def ordered_boundaries(d):
    """Innermost first. Medullary margin always innermost, outer edge (periosteum, broken edge) always outermost;
    the rest by distance from the centre, else depth below the outer edge, else kind and index (as in the tool)."""
    bs = [c for c in d["curves"] if c["kind"] in BOUNDARY and not c.get("branch_of") and segs_ok(c)]
    z = d.get("centre")
    osegs = outer_segments(d)
    if z:
        method = "centre"
        key = lambda c: sum(dist(q, z) for q in sample(c)) / len(sample(c))
    elif osegs:
        method = "periosteum_distance"
        key = lambda c: -sum(dist_to_polyline(q, osegs) for q in sample(c)) / len(sample(c))
    else:
        method = "index"
        rank = {"remodel_front": 0, "mark": 1, "label": 1, "efs_inner": 2}
        key = lambda c: rank[c["kind"]] * 1000 + (c.get("mark_index") or c.get("label_index") or 0)

    def full(c):
        if c["kind"] == "medullary":
            return (-math.inf, 0)
        if c["kind"] in OUTER:
            return (math.inf, 0 if c["kind"] == "periosteum" else 1)
        return (key(c), 0)
    return method, sorted(bs, key=full)


def zone_boundaries(order):
    """The outer edge counts once: the periosteum if traced, else the broken edge."""
    outers = [c for c in order if c["kind"] in OUTER]
    return [c for c in order if c["kind"] not in OUTER] + outers[:1]


def coverage(c_list, z):
    bins = set()
    for c in c_list:
        for s in c["segments"]:
            for q in resample(s):
                bins.add(int(angle(q, z) // BIN) % NB)
    return len(bins) * BIN


def close_periosteum(d):
    """A periosteum traced all the way round ends a click short of where it started. Close it when the gap is at most
    12% of its length and no broken edge was traced across the gap (same rule as the tool's closeLoop)."""
    brk = [c for c in d["curves"] if c["kind"] == "broken_edge"]
    out = []
    for c in d["curves"]:
        if c["kind"] == "periosteum":
            segs = []
            for s in c["segments"]:
                if len(s) >= 4:
                    a, b = s[0], s[-1]
                    gap, L = dist(a, b), sum(dist(p, q) for p, q in zip(s, s[1:]))
                    mid = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]
                    if 0 < gap <= 0.12 * L and not any(dist_to_polyline(mid, x["segments"]) < gap for x in brk):
                        s = s + [a]
                segs.append(s)
            c = dict(c, segments=segs)
        out.append(c)
    return dict(d, curves=out)


# ---------- one export ----------
def process(path, d, args):
    d = close_periosteum(d)
    issues = []
    ex = os.path.splitext(os.path.basename(path))[0]
    task = dict(d.get("task") or {})
    task.update({k: v for k, v in (d.get("identity_edits") or {}).items() if v not in (None, "")})
    k, unit = scale(d)
    if k is None:
        issues.append(f"no usable scale (unit {unit})" if unit else "no scale")
    mm = (lambda px: px / k) if k else (lambda px: None)
    z = d.get("centre")
    by_id = {c["id"]: c for c in d["curves"]}
    method, order = ordered_boundaries(d)
    rank = {c["id"]: i + 1 for i, c in enumerate(order)}

    peri = next((c for c in order if c["kind"] == "periosteum"), None)
    osegs = outer_segments(d)
    peri_bins = binned_radius(peri, z) if (peri and z) else None
    peri_cov = arc(peri_bins)[0] if peri_bins else 0
    peri_circ = star_perimeter(fill_circular(peri_bins)) if peri_bins and peri_cov >= 330 else None

    base = {"export_id": ex, "figure_id": task.get("figure_id"), "measurer": d.get("measurer"),
            "source_id": task.get("source_id"), "taxon": task.get("taxon"), "specimen": task.get("specimen"),
            "element": task.get("element")}
    curves, zones, counts = [], [], []
    for c in d["curves"]:
        if not segs_ok(c) or c["kind"] in ("ignore", "transect"):
            continue
        closed = c["kind"] == "count_area"
        L = sum(dist(a, b) for s in c["segments"] for a, b in zip(s, s[1:]))
        if closed:
            L += sum(dist(s[-1], s[0]) for s in c["segments"] if len(s) > 2)
        parent = by_id.get(c.get("branch_of")) if c.get("branch_of") else None
        if c.get("branch_of") and not parent:
            issues.append(f"branch {c['id']} points to a missing mark")
        row = dict(base, curve_id=c["id"], kind=c["kind"], radial_rank=rank.get(c["id"]),
                   mark_index=(parent or c).get("mark_index") if c["kind"] == "mark" else None,
                   is_branch=bool(parent), branch_of=c.get("branch_of"),
                   mark_type=c.get("mark_type"), confidence=c.get("confidence"), index_uncertain=c.get("index_uncertain"),
                   in_efs=c.get("in_efs"), label_index=c.get("label_index"), dye=c.get("dye"),
                   label_date=c.get("label_date") or None, age_days=c.get("age_days"),
                   n_segments=len(c["segments"]), n_points=sum(len(s) for s in c["segments"]),
                   length_mm=r4(mm(L)))
        if closed:
            counts.append((c, row))
            continue
        if z and row["n_points"] > 1:
            pts = [q for s in c["segments"] for q in resample(s)]
            rs = [dist(q, z) for q in pts]
            bins = binned_radius(c, z)
            cov, a0, a1 = arc(bins)
            row.update(radius_mean_mm=r4(mm(sum(rs) / len(rs))), radius_min_mm=r4(mm(min(rs))), radius_max_mm=r4(mm(max(rs))),
                       coverage_deg=cov, arc_start_deg=a0, arc_end_deg=a1)
            if cov >= args.min_coverage:
                row["circ_polar_mm"] = r4(mm(star_perimeter(fill_circular(bins))))
            if peri_circ and c is not peri:
                rat = [bins[i] / peri_bins[i] for i in range(NB) if bins[i] is not None and peri_bins[i]]
                if rat:
                    row["circ_peri_scaled_mm"] = r4(mm(peri_circ * sum(rat) / len(rat)))
            elif c is peri and peri_circ:
                row["circ_peri_scaled_mm"] = r4(mm(peri_circ))
        elif osegs and c["kind"] not in OUTER and row["n_points"] > 0:   # wedge with no centre: depth below the outer edge
            ds = [dist_to_polyline(q, osegs) for s in c["segments"] for q in resample(s)]
            row["depth_mean_mm"] = r4(mm(sum(ds) / len(ds)))
        curves.append(row)

    def zone_row(zd, **kw):
        ob = zd.get("obs") or {}
        r = dict(base, zone_rank=len(zones) + 1, **kw)
        r.update({f: zd.get(f) for f in DESC})
        r.update({f"obs_{f}": ob.get(f) for f in OBS})
        r.update(varies=zd.get("varies"), legacy=zd.get("legacy"),
                 described=any(zd.get(f) for f in ("matrix", "vasc_density", "vasc_orient", "remodelling", "tissue")))
        return r
    zb = zone_boundaries(order)
    for i, c in enumerate(zb):
        if c["kind"] == "medullary":
            continue
        inner = zb[i - 1] if i > 0 else None
        zones.append(zone_row(c.get("zone") or {}, inner_curve_id=inner["id"] if inner else None,
                              inner_kind=inner["kind"] if inner else None,
                              inner_mark_index=inner.get("mark_index") if inner else None,
                              outer_curve_id=c["id"], outer_kind=c["kind"], outer_mark_index=c.get("mark_index")))
    if not zones and d.get("field"):
        zones.append(zone_row(d["field"], inner_curve_id=None, inner_kind="field", inner_mark_index=None,
                              outer_curve_id=None, outer_kind="field", outer_mark_index=None))

    # transects: crossings ordered from the outer end; marks numbered 1 = innermost mark crossing on the transect
    tr_rows = []
    for ti, c in enumerate([c for c in d["curves"] if c["kind"] == "transect" and c["segments"] and len(c["segments"][0]) == 2]):
        a, b = c["segments"][0]
        xs = sorted(c.get("crossings") or [], key=lambda x: x.get("t", 0))
        n_m = sum(x.get("kind") == "mark" for x in xs)
        mi, prev = n_m + 1, None
        ends = c.get("ends") or {}
        pts = [("outer_end", a, None)] + [(x.get("kind"), x["p"], x) for x in xs] + [("inner_end", b, None)]
        rows = []
        for j, (kind, p, x) in enumerate(pts):
            idx = None
            if kind == "mark":
                mi -= 1
                idx = x.get("mark_index") if x.get("index_auto") is False else mi
            rows.append(dict(base, transect_id=c["id"], transect_n=ti + 1, order_from_outer=j, kind=kind,
                             end_type=ends.get("outer") if kind == "outer_end" else ends.get("inner") if kind == "inner_end" else None,
                             mark_index=idx, dist_from_outer_mm=r4(mm(dist(p, a))),
                             radius_mm=r4(mm(dist(p, z))) if z else None,
                             angle_deg=round(angle(p, z), 1) if z else None, x_px=round(p[0], 1), y_px=round(p[1], 1)))
        for j in range(len(rows) - 1):
            rows[j]["width_to_next_inward_mm"] = r4(rows[j + 1]["dist_from_outer_mm"] - rows[j]["dist_from_outer_mm"]) if k else None
        if not n_m:
            issues.append(f"transect {ti + 1} has no mark crossings")
        tr_rows += rows

    lac = d.get("lacunae") or []
    inside = set()
    lac_rows = []
    for c, row in counts:
        poly = c["segments"][0]
        if len(poly) < 3:
            continue
        idx = [j for j, q in enumerate(lac) if in_poly(q, poly)]
        inside.update(idx)
        a = poly_area(poly) / k / k if k else None
        lac_rows.append(dict(base, area_id=c["id"], n_lacunae=len(idx), area_mm2=r4(a),
                             density_per_mm2=r4(len(idx) / a) if a else None))
    if len(lac) - len(inside):
        issues.append(f"{len(lac) - len(inside)} lacunae outside any counting area")

    # features: recompute area, centroid and shape from the pixel runs
    feat_rows, fz_rows = [], []
    fs = d.get("features") or {}
    prof = [(c["id"], binned_radius(c, z)) for c in zb] if z else []

    def zone_of(x, y):
        if not z:
            return None
        b, r = int(angle((x, y), z) // BIN) % NB, dist((x, y), z)
        best = min(((pr[b], cid) for cid, pr in prof if pr[b] is not None and pr[b] >= r), default=None)
        return best[1] if best else None
    for o in fs.get("objects", []):
        runs = o.get("runs") or []
        n = sx = sy = sxx = syy = sxy = 0
        for i in range(0, len(runs), 3):
            y, x0, x1 = runs[i], runs[i + 1], runs[i + 2]
            m = x1 - x0 + 1
            sxs = (x0 + x1) * m / 2
            n += m; sx += sxs; sy += y * m; sxy += y * sxs; syy += y * y * m
            sxx += (x1 * (x1 + 1) * (2 * x1 + 1) - (x0 - 1) * x0 * (2 * x0 - 1)) / 6
        if not n:
            issues.append(f"feature {o.get('id')} has no pixels")
            continue
        cx, cy = sx / n, sy / n
        vxx, vyy, vxy = sxx / n - cx * cx + 1 / 12, syy / n - cy * cy + 1 / 12, sxy / n - cx * cy
        tr, det = vxx + vyy, vxx * vyy - vxy * vxy
        disc = math.sqrt(max(tr * tr / 4 - det, 0))
        l1, l2 = tr / 2 + disc, max(tr / 2 - disc, 1e-6)
        th = 0.5 * math.atan2(2 * vxy, vxx - vyy)
        zid = zone_of(cx, cy)
        feat_rows.append(dict(base, feature_id=o.get("id"), label=o.get("label"), source=o.get("source"),
                              mask_value=o.get("mask_value"), cx_px=round(cx, 1), cy_px=round(cy, 1), area_px=n,
                              area_mm2=r4(n / k / k) if k else None,
                              equiv_diameter_um=r4(2 * math.sqrt(n / math.pi) / k * 1000) if k else None,
                              elongation=r4(1 - math.sqrt(l2 / l1)), orientation_deg=round(math.degrees(th), 1),
                              radial_cos=r4(abs(math.cos(th - math.atan2(cy - z[1], cx - z[0])))) if z else None,
                              zone_outer_curve_id=zid, zone_outer_kind=by_id[zid]["kind"] if zid else None,
                              zone_outer_mark_index=by_id[zid].get("mark_index") if zid else None))
    dfe = (d.get("derived") or {}).get("features") or {}
    u2 = (dfe.get("unit2") or "")[:-1]
    f2 = TO_MM[u2] ** 2 if u2 in TO_MM else None
    for zr in dfe.get("by_zone", []):
        for lab, e in (zr.get("classes") or {}).items():
            fz_rows.append(dict(base, zone_outer_curve_id=zr.get("zone_outer_id"), zone_outer=zr.get("zone_outer"),
                                region=dfe.get("region"), label=lab, n=e.get("n"),
                                area_mm2=r4(e["area"] * f2) if f2 and e.get("area") is not None else None,
                                zone_region_area_mm2=r4(zr["region_area"] * f2) if f2 and zr.get("region_area") is not None else None,
                                area_fraction=e.get("area_fraction"),
                                density_per_mm2=r4(e["n"] / (zr["region_area"] * f2)) if f2 and zr.get("region_area") else None))

    # cross-check against the tool's own derived block (lengths only; same unit as the export)
    dv = d.get("derived")
    if dv and k and dv.get("px_per_unit"):
        f = TO_MM[dv.get("unit", "mm")] if dv.get("unit") in TO_MM else None
        mine = {r["curve_id"]: r["length_mm"] for r in curves + [r for _, r in counts]}
        for dc in dv.get("curves", []):
            v = dc.get("length", dc.get("perimeter"))
            if f and v is not None and mine.get(dc["id"]) is not None and abs(v * f - mine[dc["id"]]) > 1e-3 * max(1, abs(v * f)):
                issues.append(f"length of {dc['id']} differs from the tool's ({v * f:.4f} vs {mine[dc['id']]:.4f} mm)")

    sec = d.get("section") or {}
    marks = [r for r in curves if r["kind"] == "mark" and not r["is_branch"]]
    idxs = [r["mark_index"] for r in marks if r["mark_index"] is not None]
    if len(idxs) != len(set(idxs)):
        issues.append("duplicate mark indices")
    su = d.get("setup") or {}
    comp = su.get("companion") or {}
    cnt = d.get("counts") or {}
    tr_marks = max([sum(1 for r in tr_rows if r["transect_id"] == tid and r["kind"] == "mark")
                    for tid in {r["transect_id"] for r in tr_rows}] or [0])
    vis = cnt.get("visible_marks")
    lo = max(len(marks_main := [c for c in d["curves"] if c["kind"] == "mark" and not c.get("branch_of") and segs_ok(c)]),
             tr_marks, vis or 0)
    hi = vis if vis is not None and not cnt.get("visible_at_least") and vis >= lo else None
    if vis is not None and not cnt.get("visible_at_least") and vis < lo:
        issues.append(f"count {vis} is below the marks traced or crossed ({lo})")
    peri_c = [c for c in d["curves"] if c["kind"] == "periosteum" and segs_ok(c)]
    brk_c = [c for c in d["curves"] if c["kind"] == "broken_edge" and segs_ok(c)]
    pcov, bcov = (coverage(peri_c, z), coverage(brk_c, z)) if z else (None, None)
    steps = d.get("steps") or {}
    srow = dict(base, citation=task.get("citation"), view=su.get("view") or task.get("view"), light=su.get("light") or task.get("light"),
                task_view=task.get("view"), task_light=task.get("light"),
                companion_file=comp.get("file"), companion_light=comp.get("light"), companion_alignment=comp.get("transform_method"),
                centre_method=d.get("centre_method") or ("clicked" if z else None),
                centre_fit_method=(d.get("centre_fit") or {}).get("method"), centre_fit_arc_deg=(d.get("centre_fit") or {}).get("arc_deg"),
                n_marks_traced=len(marks_main), n_transects=len({r["transect_id"] for r in tr_rows}), max_marks_on_transect=tr_marks,
                count_visible=vis, count_at_least=cnt.get("visible_at_least"), count_lo=lo if (lo or vis is not None) else None, count_hi=hi,
                efs_lines_min=cnt.get("efs_lines_min"),
                has_broken_edge=bool(brk_c), has_remodel_front=any(c["kind"] == "remodel_front" and segs_ok(c) for c in d["curves"]),
                periosteum_cov_deg=pcov, broken_edge_cov_deg=bcov,
                # whole sections with a centre: share of the circle without periosteum (broken edges are often radial
                # cracks, so their own coverage understates the loss); otherwise broken edge vs periosteum coverage
                outer_broken_share=(r4(1 - pcov / 360) if z and (su.get("view") or "whole") == "whole" else
                                    r4(bcov / (pcov + bcov)) if z and (pcov or bcov) else
                                    1.0 if brk_c and not peri_c else 0.0 if peri_c else None),
                rules_version=d.get("rules_version"),
                steps_skipped="; ".join(f"{k}: {v.get('reason')}" for k, v in steps.items() if (v or {}).get("status") == "skipped") or None,
                quality_tier=task.get("quality_tier"), repeat=d.get("repeat"), schema_version=d.get("schema_version"),
                migrated_from=d.get("migrated_from"), tool_version=d.get("tool_version"),
                exported_at=d.get("exported_at"), time_spent_s=d.get("time_spent_s"),
                image_file=(d.get("image") or {}).get("file"), image_sha256=(d.get("image") or {}).get("sha256"),
                sha256_matches_manifest=(d.get("image") or {}).get("sha256_matches_manifest"),
                calib_unit=(d.get("calibration") or {}).get("unit"), px_per_mm=r4(k), has_centre=bool(z),
                zone_order_method=method, position=sec.get("position"), plane=sec.get("plane"),
                min_circumference=sec.get("min_circumference"), n_marks=len(marks),
                n_branches=sum(r["is_branch"] for r in curves), n_labels=sum(r["kind"] == "label" for r in curves),
                has_periosteum=peri is not None, periosteum_coverage_deg=peri_cov if peri_bins else None,
                n_zones=len(zones), n_zones_described=sum(r["described"] for r in zones), n_lacunae=len(lac),
                n_features=len(feat_rows), feature_region=dfe.get("region"),
                n_ignore_boxes=sum(c["kind"] == "ignore" for c in d["curves"]),
                based_on=(d.get("based_on") or {}).get("file"), mask_file=(fs.get("mask") or {}).get("file"))
    # flags and scores are open-ended dictionaries: flatten them with prefixes
    for kf, v in (d.get("flags") or {}).items():
        srow[f"flag_{kf}"] = v
    for kf, v in (d.get("scores") or {}).items():
        srow[f"score_{kf}"] = v
    srow["notes"] = (d.get("notes") or "").replace("\n", " ").strip() or None
    srow["validation_warnings"] = "; ".join(d.get("validation_warnings") or []) or None
    srow["import_issues"] = "; ".join(issues) or None
    return srow, curves, zones, tr_rows, lac_rows, feat_rows, fz_rows


def write(path, rows):
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="+", help="export files or folders")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "..", "data_clean"))
    ap.add_argument("--prefix", default="digitized_")
    ap.add_argument("--min-coverage", type=float, default=180,
                    help="degrees of arc a curve needs before its polar circumference is estimated (default 180)")
    args = ap.parse_args(argv)

    docs, skipped = [], []
    for p in find_exports(args.paths):
        try:
            d = json.load(open(p))
        except (OSError, json.JSONDecodeError) as e:
            skipped.append((p, f"unreadable: {e}"))
            continue
        if not isinstance(d, dict) or "curves" not in d:
            continue                                     # not an export
        d, err = migrate(d)
        if err:
            skipped.append((p, err))
            continue
        docs.append((p, d))

    latest, superseded = {}, []
    for p, d in docs:
        key = ((d.get("task") or {}).get("figure_id"), d.get("measurer"), bool(d.get("repeat")))
        if key in latest:
            old = latest[key]
            keep, drop = ((p, d), old) if (d.get("exported_at") or "") > (old[1].get("exported_at") or "") else (old, (p, d))
            latest[key] = keep
            superseded.append(drop[0])
        else:
            latest[key] = (p, d)

    S, C, Z, T, L, FT, FZ = [], [], [], [], [], [], []
    for p, d in sorted(latest.values(), key=lambda x: x[0]):
        s, c, z, t, l, ft, fz = process(p, d, args)
        S.append(s); C += c; Z += z; T += t; L += l; FT += ft; FZ += fz

    os.makedirs(args.out, exist_ok=True)
    for name, rows in (("sections", S), ("curves", C), ("zones", Z), ("transects", T), ("lacunae", L), ("features", FT),
                       ("feature_zones", FZ)):
        write(os.path.join(args.out, f"{args.prefix}{name}.csv"), rows)

    print(f"{len(S)} exports imported -> {os.path.normpath(args.out)}/{args.prefix}*.csv "
          f"({len(C)} curves, {len(Z)} zones, {len(T)} transect points, {len(L)} counting areas, {len(FT)} features)")
    for p in superseded:
        print(f"  superseded by a later export: {p}")
    for p, why in skipped:
        print(f"  skipped {p}: {why}")
    for s in S:
        if s["import_issues"]:
            print(f"  {s['export_id']}: {s['import_issues']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
