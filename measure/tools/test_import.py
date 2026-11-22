#!/usr/bin/env python3
"""Check import_exports.py against the synthetic fixtures' known geometry.

Builds exports from measure/fixtures/synthetic/truth.json by sampling the true ellipses (no browser needed):
a complete section with a centre, the same section re-exported later (must supersede the first), a partial
section traced over a 150-degree arc only, a wedge with no centre, a 0.1-schema draft with a paired double
mark, and two 0.4 exports: an eroded section (periosteum over 200 degrees, broken edge over the rest, a
remodelling front, a transect, an interval count, guided zone answers, an other-light image) and a close-up
described as one field. Imports them into a temporary folder and compares circumferences with the true ellipse perimeters.

Usage: python3 measure/tools/test_import.py
"""
import csv, json, math, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import import_exports as imp  # noqa: E402

truth = {t["figure_id"]: t for t in json.load(open(os.path.join(HERE, "..", "fixtures", "synthetic", "truth.json")))}
T = truth["syn_full_drift"]


def ell(e, a0=0.0, a1=2 * math.pi, n=48):
    out = []
    for i in range(n):
        t = a0 + (a1 - a0) * i / (n - 1 if a1 - a0 < 2 * math.pi else n)
        out.append([round(e["cx"] + e["rx"] * math.cos(t) * math.cos(e["rot"]) - e["ry"] * math.sin(t) * math.sin(e["rot"]), 1),
                    round(e["cy"] + e["rx"] * math.cos(t) * math.sin(e["rot"]) + e["ry"] * math.sin(t) * math.cos(e["rot"]), 1)])
    return out


def ramanujan(e):
    a, b = e["rx"], e["ry"]
    h = ((a - b) / (a + b)) ** 2
    return math.pi * (a + b) * (1 + 3 * h / (10 + math.sqrt(4 - 3 * h)))


def export(fid, measurer, stamp, centre=True, a0=0.0, a1=2 * math.pi, schema="0.2"):
    curves = [{"id": "peri", "kind": "periosteum", "segments": [ell(T["periosteum"], a0, a1)],
               "zone": {"matrix": "parallel_fibred", "vasc_density": "low"}},
              {"id": "med", "kind": "medullary", "segments": [ell(T["medullary"], a0, a1)]}]
    for m in T["marks"]:
        curves.append({"id": f"m{m['index']}", "kind": "mark", "segments": [ell(m, a0, a1)], "mark_index": m["index"],
                       "index_uncertain": False, "mark_type": "LAG", "confidence": 3, "in_efs": False, "branch_of": None})
    sb = T["scale_bar_px"]
    return {"schema_version": schema, "tool_version": "test", "measurer": measurer, "repeat": False,
            "task": {"figure_id": fid, "taxon": "Synthetica demonstrata", "specimen": "SYN", "element": "femur"},
            "calibration": {"p1": sb[0], "p2": sb[1], "length": 1, "unit": "mm"},
            "centre": T["centre_at_periosteum"] if centre else None, "section": {}, "flags": {}, "scores": {},
            "curves": curves, "lacunae": [], "exported_at": stamp}


def main():
    tmp = tempfile.mkdtemp(prefix="import_test_")
    src = os.path.join(tmp, "exports"); os.makedirs(src)
    docs = {"full_old": export("full", "A", "2026-01-01T00:00:00Z"),
            "full_new": export("full", "A", "2026-01-02T00:00:00Z"),
            "partial": export("partial", "A", "2026-01-01T00:00:00Z", a0=0.2, a1=0.2 + math.radians(150)),
            "wedge": export("wedge", "A", "2026-01-01T00:00:00Z", centre=False, a0=0.2, a1=1.4),
            "old": export("old", "A", "2026-01-01T00:00:00Z", schema="0.1"),
            "brokenmarks": export("brokenmarks", "A", "2026-01-01T00:00:00Z")}
    for c in docs["brokenmarks"]["curves"]:            # complete periosteum, marks preserved over 150 degrees only
        if c["kind"] == "mark":
            m = T["marks"][c["mark_index"] - 1]; c["segments"] = [ell(m, 0.2, 0.2 + math.radians(150))]
    docs["full_new"]["curves"][2]["zone"] = {"matrix": "woven", "vasc_density": "high", "vasc_orient": "laminar"}
    # a 20 x 10 px feature between marks 1 and 2 (on the +x axis); its zone is the one inside mark 2
    m1, m2 = T["marks"][0], T["marks"][1]
    x0 = int(T["centre_at_periosteum"][0] + (m1["cx"] - T["centre_at_periosteum"][0]) + (m1["rx"] + m2["rx"]) / 2 * math.cos(m1["rot"]))
    y0 = int(m1["cy"] + (m1["rx"] + m2["rx"]) / 2 * math.sin(m1["rot"]))
    runs = [v for y in range(y0, y0 + 10) for v in (y, x0, x0 + 19)]
    docs["full_new"]["features"] = {"settings": {}, "objects": [{"id": "f1", "label": "canal", "source": "manual", "runs": runs, "mask_value": 1}]}
    old = docs["old"]; old["flags"] = {"section_oblique": True}
    old["curves"].append({"id": "m3b", "kind": "mark", "segments": [ell(T["marks"][2], 1.0, 2.0)], "mark_type": "double_multiple",
                          "linked_to": "m3", "mark_index": 3})
    # ---- schema 0.4: eroded section ----
    er = export("eroded", "A", "2026-01-01T00:00:00Z", schema="0.4")
    a_cut = math.radians(200)
    er["curves"][0]["segments"] = [ell(T["periosteum"], 0, a_cut)]
    shrink = dict(T["periosteum"], rx=T["periosteum"]["rx"] * 0.97, ry=T["periosteum"]["ry"] * 0.97)
    er["curves"].append({"id": "brk", "kind": "broken_edge", "segments": [ell(shrink, a_cut, 2 * math.pi)]})
    rf = dict(T["medullary"]); m1 = T["marks"][0]
    rf.update(rx=(rf["rx"] + m1["rx"]) / 2, ry=(rf["ry"] + m1["ry"]) / 2)
    er["curves"].append({"id": "rem", "kind": "remodel_front", "segments": [ell(rf)]})
    er["curves"][0]["zone"] = {"obs": {"canal_amount": "many", "canal_pattern": "circ", "prim_osteons": "yes", "sec_osteons": "none",
                                       "matrix_xpl": "fl"}, "matrix": "woven", "vasc_density": "moderate", "vasc_orient": "laminar",
                               "remodelling": "none", "prim_osteons": True, "tissue": "fibrolamellar", "basis_light": "crossed_polars",
                               "rules_version": "obs-1"}
    for c in er["curves"]:                              # marks auto-numbered by the tool: drop the old zone fields
        c.pop("zone", None) if c["id"] != "peri" else None
    # transect along the +x direction from the periosteum to the medullary margin, crossing every mark
    zc = T["centre_at_periosteum"]
    def on_ray(e):                                       # point of ellipse e nearest the ray from zc along +x
        pts = ell(e, 0, 2 * math.pi, n=720)
        return min((p for p in pts if p[0] > zc[0]), key=lambda p: abs(p[1] - zc[1]))
    a, b = on_ray(T["periosteum"]), on_ray(T["medullary"])
    L2 = (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2
    xs = []
    for m in T["marks"]:
        q = on_ray(m); t = ((q[0] - a[0]) * (b[0] - a[0]) + (q[1] - a[1]) * (b[1] - a[1])) / L2
        xs.append({"p": [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])], "t": t, "kind": "mark", "index_auto": True})
    er["curves"].append({"id": "tr1", "kind": "transect", "segments": [[a, b]], "crossings": xs[::-1],
                         "ends": {"outer": "periosteum", "inner": "medullary"}})
    er.update(setup={"view": "whole", "light": "plain", "marks_mode": "trace",
                     "companion": {"file": "eroded_xpl.png", "light": "xpl", "transform": [1, 0, 0, 0], "transform_method": "same_size"}},
              counts={"visible_marks": 9, "visible_at_least": True, "efs_lines_min": None}, centre_method="fitted",
              centre_fit={"method": "circle_fit", "arc_deg": 360, "rms_px": 1.0},
              scores={"efs": "unclear", "efs_reason": "outer surface broken or eroded", "sec_extent": "inner"},
              steps={"more": {"status": "skipped", "reason": "not asked"}}, rules_version="obs-1")
    docs["eroded"] = er
    # ---- schema 0.4: close-up, one field ----
    cu = export("closeup", "A", "2026-01-01T00:00:00Z", centre=False, schema="0.4")
    cu["curves"] = []
    cu.update(setup={"view": "detail", "light": "xpl"}, counts={}, field={"obs": {"canal_amount": "few", "canal_pattern": "dots",
              "prim_osteons": "no", "sec_osteons": "most", "matrix_xpl": "lam"}, "tissue": "secondary_haversian", "matrix": "lamellar"})
    docs["closeup"] = cu
    for name, d in docs.items():
        json.dump(d, open(os.path.join(src, f"{name}__A__x.json"), "w"))
    out = os.path.join(tmp, "out")
    imp.main([src, "--out", out])
    rd = lambda n: list(csv.DictReader(open(os.path.join(out, f"digitized_{n}.csv"))))
    S, C, Z, TR = rd("sections"), rd("curves"), rd("zones"), rd("transects")
    fails = []
    check = lambda ok, msg: None if ok else fails.append(msg)

    k = T["px_per_mm"]
    check(sorted(s["figure_id"] for s in S) == ["brokenmarks", "closeup", "eroded", "full", "old", "partial", "wedge"], "superseded export not dropped")
    full = [c for c in C if c["figure_id"] == "full"]
    errs = []
    for m in T["marks"] + [T["periosteum"]]:
        cid = f"m{m['index']}" if "index" in m else "peri"
        row = next(c for c in full if c["curve_id"] == cid)
        true = ramanujan(m) / k
        errs.append(abs(float(row["circ_polar_mm"]) - true) / true)
    print(f"full ring, polar circumference: max relative error {max(errs):.4f}")
    check(max(errs) < 0.01, "full-ring circumference off by more than 1%")
    ranks = [int(c["radial_rank"]) for c in sorted(full, key=lambda c: (c["kind"] != "medullary", int(c["mark_index"] or 99)))]
    check(ranks == sorted(ranks), f"radial order wrong: {ranks}")

    part = [c for c in C if c["figure_id"] == "partial" and c["kind"] == "mark"]
    errs = []
    for c in part:
        m = T["marks"][int(c["mark_index"]) - 1]
        check(c["circ_polar_mm"] == "", "polar circumference given for a 150-degree arc")
        check(130 <= float(c["coverage_deg"]) <= 160, f"coverage {c['coverage_deg']} for a 150-degree arc")
        check(c["circ_peri_scaled_mm"] == "", "periosteum-scaled circumference given although the periosteum is partial")
    errs = [abs(float(c["circ_peri_scaled_mm"]) - ramanujan(T["marks"][int(c["mark_index"]) - 1]) / k)
            / (ramanujan(T["marks"][int(c["mark_index"]) - 1]) / k)
            for c in C if c["figure_id"] == "brokenmarks" and c["kind"] == "mark"]
    print(f"150-degree marks inside a complete periosteum, periosteum-scaled circumference: max relative error {max(errs):.3f}")
    check(max(errs) < 0.10, "periosteum-scaled circumference off by more than 10%")   # drift bias, see README
    wedge = [c for c in C if c["figure_id"] == "wedge" and c["kind"] == "mark"]
    check(all(c["depth_mean_mm"] for c in wedge), "wedge marks lack depth below the periosteum")
    depths = [float(c["depth_mean_mm"]) for c in sorted(wedge, key=lambda c: int(c["mark_index"]))]
    check(depths == sorted(depths, reverse=True), "wedge depths not decreasing outward")
    oldS = next(s for s in S if s["figure_id"] == "old")
    check(oldS["plane"] == "oblique" and oldS["migrated_from"] == "0.1", "0.1 migration of section_oblique failed")
    br = next(c for c in C if c["curve_id"] == "m3b")
    check(br["is_branch"] == "True" and br["mark_index"] == "3", "0.1 double mark not migrated to a branch")
    z1 = next(z for z in Z if z["figure_id"] == "full" and z["outer_curve_id"] == "m1")
    check(z1["inner_kind"] == "medullary" and z1["matrix"] == "woven", "zone inside mark 1 wrong")
    check(next(z for z in Z if z["figure_id"] == "full" and z["outer_kind"] == "periosteum")["inner_curve_id"] == "m8",
          "outermost zone not bounded by the last mark")

    FT = rd("features")
    check(len(FT) == 1 and FT[0]["zone_outer_curve_id"] == "m2", f"feature zone wrong: {[(f['label'], f['zone_outer_curve_id']) for f in FT]}")
    check(len(FT) == 1 and abs(float(FT[0]["area_mm2"]) - 200 / k / k) < 1e-4, "feature area wrong")

    # ---- 0.4 checks ----
    es = next(s for s in S if s["figure_id"] == "eroded")
    share = float(es["outer_broken_share"])
    print(f"eroded: broken share of the outer edge {share:.3f} (true {160 / 360:.3f})")
    check(abs(share - 160 / 360) < 0.04, f"broken-edge share {share}")
    check(es["count_lo"] == "9" and es["count_hi"] == "", f"interval count wrong: [{es['count_lo']}, {es['count_hi']}]")
    check(es["has_remodel_front"] == "True" and es["companion_light"] == "xpl" and es["centre_method"] == "fitted", "0.4 section fields")
    check(es["score_efs"] == "unclear" and es["score_efs_reason"].startswith("outer"), "EFS reason not imported")
    check("more: not asked" in es["steps_skipped"], "skipped step reason not imported")
    ez = [z for z in Z if z["figure_id"] == "eroded"]
    check(not any(z["outer_kind"] == "broken_edge" for z in ez), "outer edge counted twice in the zones")
    check(ez[0]["inner_kind"] == "medullary" and ez[0]["outer_kind"] == "remodel_front", f"innermost zone wrong: {ez[0]['inner_kind']}->{ez[0]['outer_kind']}")
    outer = next(z for z in ez if z["outer_kind"] == "periosteum")
    check(outer["tissue"] == "fibrolamellar" and outer["obs_canal_pattern"] == "circ" and outer["vasc_orient"] == "laminar", "zone obs not imported")
    check(len(ez) == len(T["marks"]) + 2, f"eroded zone count {len(ez)}")
    tr = [r for r in TR if r["figure_id"] == "eroded"]
    mk = [r for r in tr if r["kind"] == "mark"]
    check([int(r["mark_index"]) for r in mk] == list(range(len(T["marks"]), 0, -1)), "transect marks not numbered inner = 1")
    d_out = [float(r["dist_from_outer_mm"]) for r in tr]
    check(d_out == sorted(d_out), "transect crossings not ordered from the outer end")
    w = sum(float(r["width_to_next_inward_mm"]) for r in tr if r["width_to_next_inward_mm"])
    check(abs(w - d_out[-1]) < 1e-3, "transect widths do not add up to its length")
    errs = [abs(float(r["radius_mm"]) - math.hypot(r_x - T["centre_at_periosteum"][0], r_y - T["centre_at_periosteum"][1]) / k)
            for r in mk for r_x, r_y in [(float(r["x_px"]), float(r["y_px"]))]]
    check(max(errs) < 1e-3, "transect radii wrong")
    cz = [z for z in Z if z["figure_id"] == "closeup"]
    check(len(cz) == 1 and cz[0]["outer_kind"] == "field" and cz[0]["tissue"] == "secondary_haversian", "close-up field zone missing")

    print("\n".join(fails) if fails else "all checks passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
