#!/usr/bin/env python3
"""Import Histology Measure exports (JSON, schema 0.1 or 0.2) into flat tables for analysis.

Everything is recomputed from the raw image-pixel points and the scale bar; the `derived` block in each export
is used only as a cross-check. All lengths are written in mm (um and cm exports are converted).

Usage:
    python3 measure/tools/import_exports.py EXPORTS [EXPORTS ...] [--out data_clean] [--prefix digitized_]
EXPORTS are export files or folders (searched recursively for *.json files that carry `curves`).

Writes <out>/<prefix>{sections,curves,zones,lacunae}.csv:
  sections  one row per kept export: identity, section metadata, flags, scores, scale, import issues
  curves    one row per traced curve, innermost first (`radial_rank`), with length, radii, angular coverage and
            two circumference estimates (see README, "Import")
  zones     one row per zone between consecutive boundaries, with its descriptors
  lacunae   one row per counting area: lacunae counted, area (mm2), density (per mm2)

If the same measurer exported the same image more than once (same figure_id, measurer and repeat flag), only the
latest export is kept and the others are listed as superseded. Standard library only.
"""
import argparse, csv, json, math, os, sys
from collections import defaultdict

BOUNDARY = {"periosteum", "medullary", "mark", "label", "efs_inner"}
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
    if d.get("schema_version") == "0.2":
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


def ordered_boundaries(d):
    bs = [c for c in d["curves"] if c["kind"] in BOUNDARY and not c.get("branch_of") and segs_ok(c)]
    z = d.get("centre")
    peri = next((c for c in d["curves"] if c["kind"] == "periosteum" and any(len(s) > 1 for s in c["segments"])), None)
    if z:
        method = "centre"
        key = lambda c: sum(dist(q, z) for q in sample(c)) / len(sample(c))
    elif peri:
        method = "periosteum_distance"
        key = lambda c: math.inf if c is peri else -sum(dist_to_polyline(q, peri["segments"]) for q in sample(c)) / len(sample(c))
    else:
        method = "index"
        rank = {"medullary": 0, "mark": 1, "label": 1, "efs_inner": 2, "periosteum": 3}
        key = lambda c: rank[c["kind"]] * 1000 + (c.get("mark_index") or c.get("label_index") or 0)
    return method, sorted(bs, key=key)


# ---------- one export ----------
def process(path, d, args):
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
    peri_bins = binned_radius(peri, z) if (peri and z) else None
    peri_cov = arc(peri_bins)[0] if peri_bins else 0
    peri_circ = star_perimeter(fill_circular(peri_bins)) if peri_bins and peri_cov >= 330 else None

    base = {"export_id": ex, "figure_id": task.get("figure_id"), "measurer": d.get("measurer"),
            "source_id": task.get("source_id"), "taxon": task.get("taxon"), "specimen": task.get("specimen"),
            "element": task.get("element")}
    curves, zones, counts = [], [], []
    for c in d["curves"]:
        if not segs_ok(c):
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
        elif peri and c is not peri and row["n_points"] > 0:     # wedge with no centre: depth below the periosteum
            ds = [dist_to_polyline(q, peri["segments"]) for s in c["segments"] for q in resample(s)]
            row["depth_mean_mm"] = r4(mm(sum(ds) / len(ds)))
        curves.append(row)

    for i, c in enumerate(order):
        if c["kind"] == "medullary":
            continue
        inner = order[i - 1] if i > 0 else None
        zd = c.get("zone") or {}
        zones.append(dict(base, zone_rank=len(zones) + 1, inner_curve_id=inner["id"] if inner else None,
                          inner_kind=inner["kind"] if inner else None,
                          inner_mark_index=inner.get("mark_index") if inner else None,
                          outer_curve_id=c["id"], outer_kind=c["kind"], outer_mark_index=c.get("mark_index"),
                          matrix=zd.get("matrix"), vasc_density=zd.get("vasc_density"), vasc_orient=zd.get("vasc_orient"),
                          remodelling=zd.get("remodelling"), varies=zd.get("varies"),
                          described=any(zd.get(f) for f in ("matrix", "vasc_density", "vasc_orient", "remodelling"))))

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
    srow = dict(base, citation=task.get("citation"), view=task.get("view"), light=task.get("light"),
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
                n_zones=len(zones), n_zones_described=sum(r["described"] for r in zones), n_lacunae=len(lac))
    # flags and scores are open-ended dictionaries: flatten them with prefixes
    for kf, v in (d.get("flags") or {}).items():
        srow[f"flag_{kf}"] = v
    for kf, v in (d.get("scores") or {}).items():
        srow[f"score_{kf}"] = v
    srow["notes"] = (d.get("notes") or "").replace("\n", " ").strip() or None
    srow["validation_warnings"] = "; ".join(d.get("validation_warnings") or []) or None
    srow["import_issues"] = "; ".join(issues) or None
    return srow, curves, zones, lac_rows


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

    S, C, Z, L = [], [], [], []
    for p, d in sorted(latest.values(), key=lambda x: x[0]):
        s, c, z, l = process(p, d, args)
        S.append(s); C += c; Z += z; L += l

    os.makedirs(args.out, exist_ok=True)
    for name, rows in (("sections", S), ("curves", C), ("zones", Z), ("lacunae", L)):
        write(os.path.join(args.out, f"{args.prefix}{name}.csv"), rows)

    print(f"{len(S)} exports imported -> {os.path.normpath(args.out)}/{args.prefix}*.csv "
          f"({len(C)} curves, {len(Z)} zones, {len(L)} counting areas)")
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
