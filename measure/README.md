# Histology Measure 0.4

A browser tool for recording bone-histology section images: scale, outline, growth marks, tissue and a few
whole-bone observations. You get one JSON file per image. Images never leave your computer: the page reads a local
folder in place and uploads nothing. It is static HTML, CSS and JS, with no build step and no server.

## Quick start

The tool shows only short labels: no instructions, definitions or interpretation (a research digitizer, not a
teaching tool). This section and the glossary below are the protocol.

1. Serve the folder (`python3 -m http.server 8000 -d measure`), open http://localhost:8000 in Chrome or Edge,
   type your initials in the top bar, and **Open folder** (or **Try the demo**).
2. Work through the steps in the right-hand panel. **Next ›** or `]` moves on; the rail at the top shows what is
   done (✓), skipped (–) or not needed (faded).
   - **Image**: **View** (whole / part / close-up / unusable) and **Light**. If the folder has a photo of the same
     section in other light, choose it under **Other light** and **Align** it: click a landmark on the main image,
     the same landmark on the other image, then a second pair. **P** (or ◑) switches between the two.
   - **Scale**: **Scale bar**, click both ends, check the length (prefilled from the manifest when known).
   - **Outline**: under **Outer**, trace **Periosteum** where the real outer surface survives and **Broken edge**
     where it doesn't; under **Inner**, trace **Medullary** (the marrow-cavity margin) and **Remodel front** (inner
     limit of readable primary bone) if secondary bone has wiped out the inner marks. A periosteum traced all the way
     round closes itself. **Centre**: **Click**, **Fit** (refused when the arc is under 90° or clearly not
     circular), or **None** (wedges).
   - **Marks**: **Mark** (or `N`), click along it, `Enter` to jump a crack, `Esc` to finish. Marks are numbered from
     the inside out automatically. **EFS base** is the inner limit of the external fundamental system. For faint or
     broken sections switch to **Transect**, **New**: click the outer end, the inner end, then each mark where it
     crosses the line; **Done**. **Total marks**: your total count, traced or not; tick **≥** if some may be hidden.
   - **Tissue**: answer the picture questions for each zone (see the glossary). **= inner** copies the zone inside;
     **Copy to rest** fills the undescribed zones.
   - **Bone**: EFS, how far ringed (secondary) osteons reach, ICL, medullary bone, image quality.
   - **Export**: fix anything marked ✕, then **Export** and **Next image**.
3. If you can't do a step, pick a reason under **Skip** at the bottom of the card. If you can't answer a question,
   press **Unclear** and pick a reason. "Unclear" is useful data, so never guess.

Drag to pan, wheel to zoom, `0` fits, `[` `]` change step, `?` lists all shortcuts, ☀ shows contrast and brightness.
When you are not tracing, click a line to select it, drag its points to move them, and right-click a point to delete it.

### Glossary of screen terms

| On screen | Means |
|---|---|
| Canals: None / Few / Many / Packed | How much of the bone is canals (holes and channels): almost none; a few, scattered; many, spread through it; packed, almost spongy |
| Canal pattern: Dots / Circular / Radial / Grid / Net / Mixed | What most canals look like: round dots or short ovals; long lines running around the bone; lines running outward like spokes; lines around joined by spokes; an irregular net; a mix |
| Canal halos | Do most canals have a halo of layered bone around them with no sharp outer edge (primary osteons)? |
| Ringed osteons: None / Few / Many / Most | Round osteons with a sharp outer ring that cut across the tissue around them (secondary osteons); Most = the bone is mostly made of them |
| Matrix (XPL): Grainy / One-way / Layered / Grainy + rings | Under crossed polars the bone between the canals is dark and speckled; bright with fibres running one way; fine bright and dark layers; grainy with bright rings around the canals |
| Lacunae: Plump / Flat / Both | Asked in plain light only: the osteocyte lacunae (tiny dark specks) are plump, round and pointing every way; flat and lined up in rows; or both kinds |
| Varies | The zone's tissue changes around the section |
| EFS | Just under the outer surface, several growth lines packed tightly together with few or no canals between them |
| Ringed osteons reach: Inner / Mid / Surface | How far out from the marrow cavity secondary osteons reach |
| Inner layer (ICL) | The marrow cavity is lined by a smooth band of layered bone with a sharp line on its outer side |
| Medullary bone | A separate lining of spongy, woven-looking bone inside the marrow cavity |
| Kastschenko's line (under More) | A wavy cement line where the cortex meets the inner (endosteal) bone |
| Crushed / Cracks | The section is crushed or deformed / has cracks or infill |

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
| Wedge or fragment | View: *part*. Outline whatever edges are present; centre **Fit** (needs ~90° of periosteum) or **None**; transects work well |
| Outer cortex eroded or broken | Trace it as **Broken edge**; answer EFS **Unclear** (surface broken) |
| Inner marks remodelled away | Trace the **Remodel front**; answer how far secondary osteons reach |
| Marks too faint to trace | Switch to transects; or give only **Total marks**, with **≥** |
| Close-up of tissue | View: *close-up*. Outline and marks are skipped; answer the tissue questions for the field |
| Plain and polarized photos of one section | Pick the second one under **Other light**; **Align**; it is not measured separately |
| A mark splits in two | Trace the split part as a new mark, then edit it (`T`) and set *Branch of* |
| No usable bone | View: *unusable*, with a reason, then export |

Centre fitting: the centroid of the periosteum if it closes (≥ 330°), otherwise a least-squares circle through the
periosteum (or through all outer edges if the periosteum spans < 90°). The fit is refused when the arc is under 90°,
or when the outline is clearly not circular: RMS misfit above 5% of the radius on arcs under 240° (a flattened rib
gave 10%; a 0.7-aspect ellipse over 200° gives 3%), or above 15% on longer arcs. Fits on elliptical or drifting
sections are approximate (about 4% of the radius on the demo); the export records `centre_method` and `centre_fit`
(`rms_rel` = misfit / radius).

Closed outline: a periosteum clicked all the way round ends a click short of where it started. It is treated as
closed (for coverage, length and the centroid) when the gap is at most 12% of its length and no broken edge was
traced across the gap; the importer applies the same rule.

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
