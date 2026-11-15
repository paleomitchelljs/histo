# Histology Measure (draft 0.1)

A browser tool for measuring bone-histology section images: scale bar, section centre, growth marks
traced as curves, periosteum, medullary margin, EFS inner boundary, flags and scores. It exports one JSON
file per image.

Images never leave your computer. You open a local folder in the browser, the page reads the images in
place, and nothing is uploaded. The page itself is static HTML, CSS and JS with no build step and no server.

## Try it

- Double-click `index.html` (Chrome or Edge recommended), or serve the folder:
  `python3 -m http.server 8000 -d measure`, then open http://localhost:8000.
  Over http(s) or localhost, Chrome can write exports straight into the image folder. From a double-clicked
  file the folder is opened read-only and exports go to Downloads instead.
- Click **Demo image** for a generated section, or **Open folder** and choose `measure/fixtures/synthetic/`
  (three synthetic sections with a manifest; regenerate with `python3 measure/tools/make_fixtures.py`).

Suggested order per image: enter your measurer ID, then:
1. **Scale bar**: click both ends, then type its length (pre-filled from the manifest).
2. **Centre**, or tick "No centre".
3. **+ Periosteum** and **+ Medullary margin**.
4. **+ Growth mark** for each mark, index 1 = innermost preserved.
5. Optionally **+ EFS inner boundary**.
6. Flags, scores, **Export**, **Next**.

Tracing: click along the curve. **Enter** ends a segment, so the next click starts a new segment of the
same mark (use it across cracks and breaks; gaps are never bridged). **Esc** stops. **Backspace** removes
the last point. Drag a point to move it; right-click a point to delete it. Mouse wheel zooms; drag the
background (or hold Space) to pan; **0** fits the image.

The live statistics under each curve (length, mean radius from the centre, angular coverage) are a
preview only. All analysis is done at import from the raw points.

## Folder format

    my_batch/
      manifest.json      {"schema_version": "0.1", "batch_id": "...", "tasks": [ {...}, ... ]}
      <image files>      PNG / JPEG / WebP (browsers cannot display TIFF: convert first)
      exports/           created by the tool when it can write to the folder

Task fields: `figure_id`, `file` (required); `sha256`, `source_id`, `citation`, `taxon`, `specimen`, `element`,
`view`, `light`, `quality_tier`, `measurable`, `scale_value`, `caption`, `assignment` {measurer, repeat}.
Without a manifest, every image in the folder becomes a task and the identity fields can be typed in.

## Export file (schema 0.1)

`<figure_id>__<measurer>__<timestamp>.json` contains:
- the task row from the manifest, plus `identity_edits`;
- the image file, its sha256 and whether it matches the manifest, and its pixel size;
- measurer, `repeat`, started/exported times, `time_spent_s`;
- `calibration` {p1, p2 (pixels), length, unit} and `calibration_px_per_unit`;
- `centre` [x, y] or null, and `flags`;
- `curves[]`, each {id, kind (periosteum / medullary / mark / efs_inner), segments [[[x, y], ...], ...]}. Marks also
  carry mark_index, index_uncertain, mark_type (LAG / annulus / double_multiple / uncertain), confidence 1-3,
  linked_to (id of the paired mark of a double LAG) and in_efs;
- `scores`, `notes`, `validation_warnings`.

All coordinates are image pixels, origin at the top left.

Work in progress is autosaved in the browser for convenience; the exported file is the record.

## Not yet in this draft

Tissue-region polygons, osteon and vascular-field tools, fixed-transect quick mode, edge-snapping, a JSON
Schema file and its validation, GitHub Pages deployment.
