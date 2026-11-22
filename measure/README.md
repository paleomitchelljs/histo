# Histology Measure 0.4

A browser tool for recording bone-histology section images: scale, outline, growth marks, tissue and a few
whole-bone observations. You get one JSON file per image. Images never leave your computer: the page reads a local
folder in place and uploads nothing. It is static HTML, CSS and JS, with no build step and no server.

## Quick start

The tool itself shows no instructions or interpretation (a research digitizer, not a teaching tool); this
section is the protocol.

1. Serve the folder (`python3 -m http.server 8000 -d measure`), open http://localhost:8000 in Chrome or Edge,
   type your initials in the top bar, and **Open folder** (or **Try the demo**).
2. Work through the steps in the right-hand panel. **Next ›** or `]` moves on; the rail at the top shows what is
   done (✓), skipped (–) or not needed (faded).
   - **Image**: what it shows (whole section / part / close-up / nothing usable) and how it was lit. If the folder
     has a crossed-polars photo of the same section, pick it here and **Align** it with two landmarks. **P** then
     switches between the two.
   - **Scale**: click both ends of the scale bar and check the length.
   - **Outline**: trace the *periosteum* where the real outer surface survives and the *broken edge* where it
     doesn't; trace the *medullary margin*, and the *remodelling front* if secondary bone has wiped out the inner
     marks. Then set the centre: **Click it**, **Fit from outline**, or **No centre** (wedges).
   - **Growth marks**: trace each mark (**Growth mark** or `N`, click along it, `Enter` to jump a crack, `Esc` to finish). Marks are
     numbered from the inside out automatically. For faint or broken sections, switch to **Transect**. Then give your total count, ticking "at least" if some may be hidden.
   - **Tissue**: answer the picture questions for each zone. **Copy to undescribed zones** saves repetition.
   - **Whole bone**: EFS, how far secondary osteons reach, inner circumferential layer, medullary bone, quality.
   - **Export**: fix anything marked ✕, then **Export** and **Next image**.
3. If you can't do a step, pick a reason under **Not done** at the bottom of the card. If you can't
   answer a question, press **Can't tell** and say why. "Can't tell" is useful data, so never guess.

Drag to pan, wheel to zoom, `0` fits, `[` `]` change step, `?` lists all shortcuts. When you are not tracing, click
a line to select it, drag its points to move them, and right-click a point to delete it.

## What gets recorded, and why

The export holds what the growth model needs, including *why* something is missing:

- **Outer vs inner loss of marks** (model 38): the periosteum vs broken edge (where the outer surface is gone),
  the remodelling front and secondary-osteon extent (where inner marks are wiped out), and EFS lines that can't
  be separated (`counts.efs_lines_min`).
- **EFS three-state**: present / absent / unclear, with a reason, so "surface broken" never becomes "absent".
- **Interval counts**: traced marks, transect crossings and the student's count become `[count_lo, count_hi]` at
  import (open-ended if "at least").
- **Radii and zone widths**: traced curves (with a clicked or fitted centre) and transects, both in image pixels
  plus the scale.
- **Tissue per zone**: the raw answers (`zone.obs`) and the descriptors derived from them by `guide.js`
  (`matrix`, `vasc_density`, `vasc_orient`, `remodelling`, `prim_osteons`, `tissue`, `basis_light`,
  `rules_version`). The derived terms are not shown on screen, so they cannot steer the answers. The rules can be
  revised and re-applied without measuring again. A matrix call from plain
  light (lacuna shape) is marked `basis_light = plain_light_lacunae`.

## Odd cases

| Case | What to do |
|---|---|
| Wedge or fragment | Image: *part of a section*. Outline whatever edges are present; centre **Fit** (needs ~90° of periosteum) or **No centre**; transects work well |
| Outer cortex eroded or broken | Trace it as **Broken edge**; answer EFS **Can't tell** (surface broken) |
| Inner marks remodelled away | Trace the **Remodelling front**; answer how far secondary osteons reach |
| Marks too faint to trace | Switch to transects; or give only the count, with "at least" |
| Close-up of tissue | Image: *close-up*. Outline and marks are skipped; answer the tissue questions for the field |
| Plain and polarized photos of one section | Pick the second one as the other-light image; **Align**; it is not measured separately |
| A mark splits in two | Trace the split part as a new mark, then edit it (`T`) and set *Splits off* |
| No usable bone | Image: *nothing usable*, with a reason, then export |

Centre fitting: the centroid of the periosteum if it closes (≥ 330°), otherwise a least-squares circle through the
periosteum (or through all outer edges if the periosteum spans < 90°). Fits on elliptical or drifting sections
are approximate (about 4% of the radius on the demo); the export records `centre_method` and `centre_fit`.

## More (optional) tools

Under **More**: the feature finder (canals and osteons), osteocyte lacuna counts, and ignore boxes.

### Feature finder

**Detect** marks every connected region that is darker (or brighter, or either) than the tissue around it.
Contrast is measured against the tissue within **Window px** and scaled by its local spread, so the section edge
and uneven staining matter less. **Sensitivity** sets the cut-off (default 1 local SD). Regions smaller than
**min area** are dropped. A connected network (reticular canals) larger than one tile (4 windows, at least 24
px) is cut into tiles, so each piece can be labelled. The search region is the traced Areas if there are any;
otherwise the cortex inside the periosteum (minus the medullary cavity) when both are traced most of the way
round the centre; otherwise all tissue. Ignore boxes are always left out. The slide background is the border
colour (white in plain light, black under crossed polars).

**Cartoon** (V) shows the candidates beside the image in the same view: grey = unlabelled, colours = labels,
pale colours = suggestions. Click a feature (in the cartoon, or on the image with **Label**, L) to label it, or
hover and press 1–6 (X = not a feature, U = clear). Where detection missed something, clicking empty tissue adds
a point feature. **Find similar** (F) compares each unlabelled candidate with the labelled ones on size, shape,
tone, colour, orientation relative to the centre, and the ring of tissue around it (which separates an osteon's
canal from a plain canal), and suggests the nearest label within **Match**. **Accept** (A) takes all
suggestions. Label some "not a feature" examples too: without negatives, nearly everything gets a suggestion.
Labelling is not covered by Undo; relabel instead.

The export stores the labelled features under `features.objects` (label, source: manual / similar / point,
pixel runs `[y, x0, x1, ...]`, centroid, area, shape) and the settings used. `<export>.mask.png` holds the whole
segmentation: each feature has its own value R + 256·G + 65536·B (0 = none); values listed in `objects[].mask_value`
are labelled, the rest are unlabelled candidates. It opens in Fiji or Python. **Load export** (choose the .json
and its .mask.png) reopens an export with exactly that segmentation, so a colleague can relabel features and
export again; the new file records `based_on`. `derived.features` gives counts, areas and area fractions per
label, overall and per zone (zone = nearest boundary outward from the centre).

Tested on the 18 Garros 2025 *Troodon* figures (plain light, about 8 µm per pixel). At that resolution the
canal network is picked up well; osteons are found through their canals, and their outlines are not traced.
Crossed-polar images are untested.


## Folder format

    my_batch/
      manifest.json      {"schema_version": "0.1", "batch_id": "...", "tasks": [ {...}, ... ]}
      <image files>      PNG / JPEG / WebP (convert TIFF first)
      exports/           created by the tool when it can write to the folder

Task fields: `figure_id`, `file` (required); `sha256`, `source_id`, `citation`, `taxon`, `specimen`, `element`,
`view`, `light`, `quality_tier`, `measurable`, `scale_value`, `caption`, `assignment` {measurer, repeat}.
Without a manifest, each image becomes a task and the identity fields can be typed in.

Optional task fields for paired images: `companion_file`, `companion_light`.

## Export file (schema 0.4)

`<figure_id>__<measurer>__<timestamp>.json`. All coordinates are image pixels, origin top left.

- `task` (manifest row), `identity_edits`, `image` {file, sha256, sha256_matches_manifest, px size}, `measurer`,
  `repeat`, `started_at`, `exported_at`, `time_spent_s`, `rules_version`;
- `setup` {view: whole | part | detail | unusable, light: plain | xpl | lambda | fluorescence | ct | other,
  companion {file, light, sha256, transform [a, b, tx, ty] (main = [a −b; b a]·companion + t), transform_method},
  marks_mode};
- `calibration` {p1, p2, length, unit}, `centre` [x, y] or null, `centre_method` clicked | fitted | null,
  `centre_fit` {method, arc_deg, rms_px};
- `section` {position, plane, min_circumference}, `flags`, `scores` (each question's answer plus
  `<key>_reason` for "unclear"), `counts` {visible_marks, visible_at_least, efs_lines_min}, `steps`
  {id: {status: skipped, reason}}, `notes`, `validation_warnings`;
- `curves[]`: {id, kind, segments [[[x, y], ...], ...]}. Kinds: periosteum, broken_edge, medullary,
  remodel_front, mark, label, efs_inner, transect, count_area, ignore.
  - marks: mark_index, index_auto, index_uncertain, mark_type (LAG / annulus / uncertain), confidence 1–3,
    in_efs, branch_of;
  - labels: label_index, dye, label_date, age_days;
  - transects: one segment [outer end, inner end], `ends` {outer, inner}, `crossings[]` {p, t (0 at the outer
    end), kind, mark_index, index_auto}. Marks are numbered 1 = innermost crossing on that transect;
  - boundary curves may carry `zone` (the tissue just inside): `obs` plus the derived descriptors above, `varies`,
    and `legacy` for zones described in 0.3 and not yet re-answered;
- `field`: the same block for a close-up (no zones);
- `lacunae[]`, `features` (see the feature finder), `based_on`;
- `derived` (convenience; recompute from the raw points): per-curve lengths, radii and coverage, transect
  distances and radii, zones in radial order, outer-edge coverage, step states, lacuna counts, features.

Zone order: the medullary margin is always innermost and the outer edge (periosteum, else broken edge) always
outermost; the rest go by distance from the centre, else depth below the outer edge, else kind and index.

Exports and drafts from 0.1–0.3 open and migrate: their mark numbers are kept as typed (`index_auto` false), and
their zone descriptors are kept (`legacy`) until someone answers the new questions.

## Import

    python3 measure/tools/import_exports.py <exports folder or files> [--out data_clean] [--min-coverage 180]

Writes `digitized_sections.csv`, `digitized_curves.csv`, `digitized_zones.csv`, `digitized_transects.csv`,
`digitized_lacunae.csv`, `digitized_features.csv` and `digitized_feature_zones.csv`. Everything is recomputed from
the raw points (the `derived` block is only cross-checked), lengths are in mm, and older exports are migrated. If one
measurer exported an image more than once, the latest export is kept. Problems go in `import_issues`.

New in 0.4: sections carry `view`, `light`, the other-light image, `centre_method`, `count_lo` / `count_hi`,
`efs_lines_min`, `periosteum_cov_deg`, `broken_edge_cov_deg`, `outer_broken_share` (whole sections with a centre:
share of the circle without periosteum), `has_remodel_front`, `score_*_reason` and `steps_skipped`. Zones carry
`obs_*` and the derived descriptors. `digitized_transects.csv` has one row per transect point (outer end,
crossings, inner end) with distance from the outer end, radius from the centre and width to the next point inward.

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

`python3 measure/tools/test_import.py` checks the import against the synthetic fixtures (including a 0.4 eroded
section with a transect, and a close-up).

## Not yet

A picture glossary with example crops for each answer; transect crossings that take their numbers from traced marks;
an ellipse fit for the centre; edge snapping; a JSON Schema file; a drift-aware circumference for partial arcs;
osteon outlines (not just canals); and the feature finder tested on crossed-polar images.
