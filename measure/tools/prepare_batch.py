#!/usr/bin/env python3
"""Turn a folder of figure images (TIFF etc.) into a batch the measuring tool can open.

Browsers cannot display TIFF, so each image is written as a lossless PNG (8-bit RGB). A manifest.json lists one task
per image with its sha256; shared fields (citation, taxon, scale value, light) can be given on the command line and
edited in the manifest afterwards. Needs Pillow.

Usage:
    python3 measure/tools/prepare_batch.py SRC_DIR OUT_DIR [--batch-id ID] [--citation "..."] [--taxon "..."]
           [--scale-value "1000 um"] [--light normal|polarized] [--source-id ID]
"""
import argparse, hashlib, json, os, sys

from PIL import Image

EXT = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".webp", ".bmp")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src"); ap.add_argument("out")
    for k in ("batch-id", "citation", "taxon", "scale-value", "light", "source-id"):
        ap.add_argument("--" + k)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    tasks = []
    for name in sorted(os.listdir(a.src)):
        if not name.lower().endswith(EXT):
            continue
        stem = os.path.splitext(name)[0]
        im = Image.open(os.path.join(a.src, name))
        if im.mode != "RGB":
            im = im.convert("RGB")
        dst = os.path.join(a.out, stem + ".png")
        im.save(dst, optimize=False)
        sha = hashlib.sha256(open(dst, "rb").read()).hexdigest()
        t = {"figure_id": stem, "file": stem + ".png", "sha256": sha, "original_file": name}
        for k in ("citation", "taxon", "scale_value", "light", "source_id"):
            v = getattr(a, k)
            if v:
                t[k] = v
        tasks.append(t)
        print(f"{name} -> {stem}.png ({im.width} x {im.height})")
    man = {"schema_version": "0.1", "batch_id": a.batch_id or os.path.basename(os.path.normpath(a.out)), "tasks": tasks}
    json.dump(man, open(os.path.join(a.out, "manifest.json"), "w"), indent=1)
    print(f"{len(tasks)} images, manifest at {os.path.join(a.out, 'manifest.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
