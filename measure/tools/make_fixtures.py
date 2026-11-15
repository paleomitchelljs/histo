#!/usr/bin/env python3
"""Synthetic section images with KNOWN geometry, for trying the measuring tool and testing the import.

Writes measure/fixtures/synthetic/: PNG images, manifest.json + manifest.csv (the folder format the tool
reads), and truth.json (true centre, ellipse of every growth mark, periosteum, medullary margin, EFS inner
boundary, scale in px per mm). No published figures are used.

Usage: python3 measure/tools/make_fixtures.py
"""
import csv, hashlib, json, math, os, random
import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fixtures", "synthetic")
random.seed(7); np.random.seed(7)


def ell_pts(cx, cy, rx, ry, rot, n=720):
    t = np.linspace(0, 2 * math.pi, n, endpoint=False)
    x = cx + rx * np.cos(t) * math.cos(rot) - ry * np.sin(t) * math.sin(rot)
    y = cy + rx * np.cos(t) * math.sin(rot) + ry * np.sin(t) * math.cos(rot)
    return list(zip(x.tolist(), y.tolist()))


def section(name, W, H, px_per_mm, centre, outer, medul, n_marks, rot, drift, efs_from=None, broken=None,
            contrast=1.0, blur=0.0, noise=6):
    """Concentric ellipses whose centres drift with growth. Fractions follow a saturating curve."""
    img = Image.new("L", (W, H), 232); d = ImageDraw.Draw(img)

    def shape(f):  # f = 0 medullary margin ... 1 periosteum
        return dict(cx=centre[0] + drift[0] * (1 - f), cy=centre[1] + drift[1] * (1 - f),
                    rx=medul[0] + (outer[0] - medul[0]) * f, ry=medul[1] + (outer[1] - medul[1]) * f, rot=rot)
    per = shape(1.0)
    d.polygon(ell_pts(**per), fill=int(232 - 70 * contrast))
    # vascular canals
    for _ in range(int(W * H / 400)):
        f = random.uniform(0.03, 0.98); a = random.uniform(0, 2 * math.pi); e = shape(f)
        x = e["cx"] + e["rx"] * math.cos(a) * math.cos(rot) - e["ry"] * math.sin(a) * math.sin(rot)
        y = e["cy"] + e["rx"] * math.cos(a) * math.sin(rot) + e["ry"] * math.sin(a) * math.cos(rot)
        r = random.uniform(1.2, 3.2)
        d.ellipse([x - r, y - r * 0.7, x + r, y + r * 0.7], fill=int(232 - 150 * contrast))
    # growth marks: saturating spacing (wide inner zones, narrow outer zones), EFS = tight outer marks
    t = np.arange(1, n_marks + 1)
    fr = 1 - np.exp(-0.30 * (t + 1.0)); fr = fr / fr[-1] * 0.985
    marks = []
    for i, f in enumerate(fr):
        e = shape(float(f)); d.line(ell_pts(**e) + [ell_pts(**e)[0]], fill=int(232 - 175 * contrast), width=2)
        marks.append(dict(index=i + 1, fraction=float(f), **e))
    med = shape(0.0); d.polygon(ell_pts(**med), fill=int(232 - 185 * contrast))
    efs = None
    if efs_from is not None:
        efs = dict(inner_boundary_mark=efs_from, **{k: marks[efs_from - 1][k] for k in ("cx", "cy", "rx", "ry", "rot")})
    if broken:  # remove a sector (fossil breakage)
        a0, a1 = broken; poly = [tuple(centre)]
        for a in np.linspace(a0, a1, 40):
            rr = max(W, H) * (1 + random.uniform(-0.03, 0.03)); poly.append((centre[0] + rr * math.cos(a), centre[1] + rr * math.sin(a)))
        d.polygon(poly, fill=232)
    arr = np.asarray(img, dtype=float) + np.random.normal(0, noise, (H, W))
    img = Image.fromarray(np.clip(arr, 0, 255).astype("uint8"))
    if blur: img = img.filter(ImageFilter.GaussianBlur(blur))
    d = ImageDraw.Draw(img)                                   # scale bar: 1 mm
    x0, y0 = W - 60 - px_per_mm, H - 50
    d.rectangle([x0, y0, x0 + px_per_mm, y0 + 7], fill=10)
    try: font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 22)
    except Exception: font = ImageFont.load_default()
    d.text((x0 + px_per_mm / 2 - 22, y0 - 30), "1 mm", fill=10, font=font)
    path = os.path.join(OUT, name + ".png"); img.convert("RGB").save(path, optimize=True)
    truth = dict(figure_id=name, px_per_mm=px_per_mm, scale_bar_px=[[x0, y0 + 3.5], [x0 + px_per_mm, y0 + 3.5]],
                 centre_at_periosteum=[per["cx"], per["cy"]], periosteum=per, medullary=med, marks=marks, efs=efs,
                 broken_sector_rad=broken, drift_px=drift,
                 note="ellipse parameters in image pixels; mark centres move by drift*(1 - fraction)")
    return path, truth


def main():
    os.makedirs(OUT, exist_ok=True)
    specs = [
        dict(name="syn_full_drift", W=1400, H=1100, px_per_mm=180, centre=(700, 560), outer=(460, 390), medul=(150, 120),
             n_marks=8, rot=0.35, drift=(-70, 30), view="full_section", caption="Synthetic full section: 8 marks, centre drifts 76 px with growth."),
        dict(name="syn_broken_efs", W=1300, H=1000, px_per_mm=150, centre=(640, 520), outer=(430, 360), medul=(140, 115),
             n_marks=11, rot=-0.2, drift=(40, -55), efs_from=9, broken=(-0.7, 0.9), view="partial_section",
             caption="Synthetic broken section: 11 marks, the outer 3 packed (EFS from mark 9), one sector missing."),
        dict(name="syn_faint", W=1200, H=950, px_per_mm=220, centre=(600, 480), outer=(380, 330), medul=(130, 100),
             n_marks=6, rot=0.9, drift=(-25, -20), view="full_section", contrast=0.45, blur=1.2, noise=10,
             caption="Synthetic low-contrast section: 6 faint marks (try the contrast slider)."),
    ]
    tasks, truths = [], []
    for s in specs:
        view, cap = s.pop("view"), s.pop("caption")
        path, truth = section(**s)
        sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
        tasks.append(dict(figure_id=s["name"], file=os.path.basename(path), sha256=sha, source_id="synthetic",
                          citation="Synthetic fixture (measure/tools/make_fixtures.py)", taxon="Synthetica demonstrata",
                          specimen="SYN-" + s["name"].split("_", 1)[1].upper(), element="femur", view=view, light="normal",
                          quality_tier="A" if view == "full_section" else "B", measurable="marks,medullary,efs",
                          scale_value="1 mm", caption=cap, assignment=dict(measurer=None, repeat=False)))
        truths.append(truth)
    json.dump(dict(schema_version="0.1", batch_id="synthetic_fixtures", tasks=tasks), open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    with open(os.path.join(OUT, "manifest.csv"), "w", newline="") as f:
        cols = [k for k in tasks[0] if k != "assignment"]; w = csv.DictWriter(f, cols, extrasaction="ignore"); w.writeheader(); w.writerows(tasks)
    json.dump(truths, open(os.path.join(OUT, "truth.json"), "w"), indent=1)
    print(f"wrote {len(tasks)} images + manifest.json, manifest.csv, truth.json to {OUT}")


if __name__ == "__main__":
    main()
