# Histology Measure 0.2

A browser tool for tracing bone-histology section images. You trace curves (growth marks, fluorochrome labels,
periosteum, medullary margin, EFS inner boundary), describe the tissue zones between them, and optionally count
osteocyte lacunae. It exports one JSON file per image.

Images never leave your computer. The page reads a local folder in place and uploads nothing. It is static HTML,
CSS and JS: no build step, no server.

The screen shows only what you need to act on. Lengths, radii, angular coverage, zone order and lacuna density
are computed at export and written to the file under `derived`.

## Run it

- Serve the folder: `python3 -m http.server 8000 -d measure`, then open http://localhost:8000 in Chrome or Edge.
  Exports then go straight into `<image folder>/exports/`. Double-clicking `index.html` also works, but exports
  then go to Downloads.
- **Demo** generates a practice section. **Open folder** on `measure/fixtures/synthetic/` gives three synthetic
  sections with a manifest (regenerate them with `python3 measure/tools/make_fixtures.py`).

## Per image

1. Measurer ID (top bar; remembered).
2. **Scale**: click both ends of the bar. The length is pre-filled from the manifest's `scale_value`; check it.
3. **Centre**: click it, or toggle **No centre** for wedges and fragments.
4. **Section**: position along the shaft, plane, and whether it is the minimum circumference.
5. **New** (N) for each curve, then click along it. Each new curve starts as the next growth mark (1 = innermost
   preserved). Use **Tag** (T) to change what it is:
   - *growth mark*: index, type (LAG / annulus / uncertain), confidence 1–3, inside EFS;
   - *branch*: a mark that splits off another over part of the cortex. Trace only the split part and set
     **Branch of** to the main mark. It takes the main mark's index and is drawn dashed, labelled e.g. `2b`;
   - *fluorochrome label*: index, dye, injection date and/or the animal's age in days at injection;
   - periosteum, medullary margin, EFS inner boundary.
6. **Zones** fill in once boundaries are traced, listed outer to inner (`2 → 3` is the tissue between marks 2
   and 3). Click one to set matrix, vascular density, canal orientation and secondary osteons. **Same as inner**
   copies the zone just inside. Orange dashes = not described yet.
7. **Lacunae** (optional, high-resolution images only): **Area** traces the counting polygon, **Count** places a
   point per lacuna. Right-click a point to remove it.
8. **Whole section** scores, notes, **Export**, **Next**.

Tracing: **Enter** leaves a gap so the next click starts a new segment of the same curve (use it across cracks;
gaps are never bridged). **Backspace** removes the last point, **Esc** stops. In **Edit** mode you can drag points
and the scale/centre handles; right-click deletes a point. Wheel or +/− zooms, arrows or Space-drag pan, 0 fits.
**?** lists the shortcuts.

The export check line under the notes says what is *Needed* (blocks export) and what to *Check* (exported as
warnings).

## Folder format

    my_batch/
      manifest.json      {"schema_version": "0.1", "batch_id": "...", "tasks": [ {...}, ... ]}
      <image files>      PNG / JPEG / WebP (convert TIFF first)
      exports/           created by the tool when it can write to the folder

Task fields: `figure_id`, `file` (required); `sha256`, `source_id`, `citation`, `taxon`, `specimen`, `element`,
`view`, `light`, `quality_tier`, `measurable`, `scale_value`, `caption`, `assignment` {measurer, repeat}.
Without a manifest, each image becomes a task and the identity fields can be typed in.

## Export file (schema 0.2)

`<figure_id>__<measurer>__<timestamp>.json`:

- `task` (manifest row), `identity_edits`, `image` {file, sha256, sha256_matches_manifest, px size};
- `measurer`, `repeat`, `started_at`, `exported_at`, `time_spent_s`;
- `calibration` {p1, p2 in pixels, length, unit}, `centre` [x, y] or null;
- `section` {position, plane, min_circumference}, `flags`, `scores`, `notes`, `validation_warnings`;
- `curves[]`: {id, kind, segments [[[x, y], ...], ...]} with kind one of periosteum, medullary, mark, label,
  efs_inner, count_area.
  - marks: mark_index, index_uncertain, mark_type, confidence, in_efs, branch_of (id of the main mark or null);
  - labels: label_index, dye, label_date, age_days;
  - any boundary curve may carry `zone` {matrix, vasc_density, vasc_orient, remodelling, varies}: the tissue
    just *inside* that curve;
- `lacunae[]`: [x, y] points;
- `derived` (convenience; recompute from the raw points for analysis):
  - `unit`, `px_per_unit`;
  - per curve: length (perimeter for count areas), radius mean/min/max from the centre, `coverage_deg` (5° bins)
    and `arc_deg` [start, end], clockwise from the +x axis in image coordinates;
  - `zones[]` in radial order with inner/outer curve ids and descriptors. `reference_for_zone_order` says how
    the order was found: distance from the centre, else from the periosteum (wedges), else kind and index;
  - `lacuna_counts[]` per area: n, area (unit²), density; `lacunae_outside_areas`.

All coordinates are image pixels, origin top left. Work in progress autosaves in the browser; the exported file
is the record. Drafts saved by 0.1 are migrated on open: paired double marks become branches, and "section
oblique" becomes plane = oblique.

## Import

    python3 measure/tools/import_exports.py <exports folder or files> [--out data_clean] [--min-coverage 180]

Writes `digitized_sections.csv`, `digitized_curves.csv`, `digitized_zones.csv` and `digitized_lacunae.csv`.
Everything is recomputed from the raw points (the `derived` block is only cross-checked), lengths are in mm,
and 0.1 exports are migrated. If one measurer exported an image more than once, the latest export is kept.
Problems (no scale, duplicate mark indices, lacunae outside a counting area, disagreement with the tool's own
lengths) go in `import_issues` and are printed.

Each curve has two circumference estimates, and both need a centre:

- `circ_polar_mm`: the radius is binned in 5° steps from the centre, gaps are interpolated around the circle,
  and the perimeter is measured. It is given only when the curve covers at least `--min-coverage` degrees.
  Error is 0.2% on complete synthetic rings.
- `circ_peri_scaled_mm`: the periosteum's perimeter multiplied by the mean ratio of the curve's radius to the
  periosteum's radius over the bins the curve covers. It is given only when the periosteum covers at least
  330°. It is meant for partly preserved inner marks. It is biased upward by centre drift, because inner
  marks are off-centre relative to the periosteum: on the drifting fixture (76 px of drift on a 460 px
  radius), a mark traced over 150° comes out 0.1–7.6% too large, worst for the innermost mark. Ratios of
  arc length do no better. Treat partial-arc circumferences as noisier, using `coverage_deg`.

Wedges with no centre get `depth_mean_mm`, the mean distance below the periosteum.

`python3 measure/tools/test_import.py` builds exports from the synthetic fixtures' true ellipses and checks the
import against them.

## Not yet

Edge snapping, fixed-transect quick mode, a JSON Schema file, a drift-aware circumference for partial arcs,
and a way to describe how a zone varies around the section other than the notes.
