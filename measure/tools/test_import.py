#!/usr/bin/env python3
"""Check import_exports.py against the synthetic fixtures' known geometry.

Builds exports from measure/fixtures/synthetic/truth.json by sampling the true ellipses (no browser needed):
a complete section with a centre, the same section re-exported later (must supersede the first), a partial
section traced over a 150-degree arc only, a wedge with no centre, and a 0.1-schema draft with a paired double
mark. Imports them into a temporary folder and compares circumferences with the true ellipse perimeters.

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
    old = docs["old"]; old["flags"] = {"section_oblique": True}
    old["curves"].append({"id": "m3b", "kind": "mark", "segments": [ell(T["marks"][2], 1.0, 2.0)], "mark_type": "double_multiple",
                          "linked_to": "m3", "mark_index": 3})
    for name, d in docs.items():
        json.dump(d, open(os.path.join(src, f"{name}__A__x.json"), "w"))
    out = os.path.join(tmp, "out")
    imp.main([src, "--out", out])
    rd = lambda n: list(csv.DictReader(open(os.path.join(out, f"digitized_{n}.csv"))))
    S, C, Z = rd("sections"), rd("curves"), rd("zones")
    fails = []
    check = lambda ok, msg: None if ok else fails.append(msg)

    k = T["px_per_mm"]
    check(sorted(s["figure_id"] for s in S) == ["brokenmarks", "full", "old", "partial", "wedge"], "superseded export not dropped")
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

    print("\n".join(fails) if fails else "all checks passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
