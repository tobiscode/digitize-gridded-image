# MEMORY.md — Working Notes for AI Agents

This file is the living scratch pad for AI coding agents working on this
repository. It records the current goal, the active plan, and verified facts
about the figures and the environment. Keep it up to date: add new findings,
remove statements that turned out to be wrong or obsolete, and update the
plan status as work progresses.

## Current goal

Add automatic handling of multiple images (branch `auto-multi-image`): the two
figures `arney_fig26.png` and `arney_fig27.png` (from Arney et al. 2014) each
contain a grid of small gridded-image panels. Each row of panels shares one
colorbar at the far right, and each column of panels shares one data grid.
Manually reading off grid lengths, corners, and colorbar tick positions is
impractical because there are no gridlines and masked (missing) data is simply
blank. We need a script that, given an approximate bounding box of the panel
array and the number of rows and columns, determines the best-fitting exact
grid (shared by all panels) and the colorbar extents, and then extracts the
data values for every panel using the logic of
`digitize_gridded_image.py`.

## Design decisions (confirmed with the user)

- Input geometry: one approximate bounding box of the whole panel array plus
  `n_rows` and `n_cols`.
- Output: full data extraction for every panel *and* the detected geometry
  (saved as JSON, plus an optional debug overlay image).
- Colorbar values: the current figures have *no tick marks*, so there is
  nothing to detect. Instead, the user provides the min/max values of each
  colorbar (one pair per row) separately; the script maps the two ends of the
  colorbar to those values. Which end holds min vs max must be explicit in the
  input (e.g. values given in top-to-bottom order) because the printed
  orientation is not readable from pixels; the debug output shows the color at
  each colorbar end so the user can sanity-check the orientation.
- Panels are delimited only by white background gaps (no axis frames).
- Masked cells (blank/fill-colored cells) are reported as NaN.
- Fully shared geometry: rows are consistent across columns and columns
  across rows, so every panel has the same shape and every cell grid is
  aligned with all others. The layout is therefore one regular lattice, and
  exactly ONE cell grid (nx, ny and cell size) is detected jointly such that
  it best matches all cells of all panels at once - not per-panel detection.
- One colorbar at the far right of each row.

## Plan

Status markers: [ ] not started, [~] in progress, [x] done.

1. [ ] Refactor `digitize_gridded_image.py` into importable functions:
   move the `__main__` pipeline into functions (image loading, per-cell
   dominant-color extraction, colorbar extraction, color-to-value mapping,
   saving, image recreation). Keep `geometric_median` and
   `project_point_on_line_3d` untouched, and keep the command-line interface
   (argument names and semantics) exactly as documented in `README.md`.
2. [ ] Numerically verify the structure of `arney_fig26.png` /
   `arney_fig27.png` (dominant colors, white-gap row and column runs,
   colorbar location) to calibrate the detection parameters.
3. [ ] New module `grid_detection.py`:
   - `detect_background_color`: modal near-white color inside the approximate
     box, doubling as the mask-fill color detector.
   - `detect_panel_lattice`: because rows are consistent across columns (and
     vice versa), the whole array is one regular lattice. Detect it as a whole
     (background-gap lines, then a least-squares/cell-size fit so that all
     panels get identical shapes and aligned column and row boundaries),
     rather than detecting rows and columns independently.
   - `detect_cell_grid`: ONE shared cell grid for all panels: score candidate
     cell counts (nx, ny) by the color-discontinuity evidence along the
     internal cell boundaries of a uniform lattice, evaluated jointly over
     all panels and aligned with the shared panel lattice (a single (nx, ny)
     and cell size that best matches every cell everywhere); optional manual
     override.
   - `detect_colorbar_ends`: locate the thin vertical color strip in the
     right margin of each row band and determine the pixel positions of its
     top and bottom ends (no ticks exist in these figures); the data values
     at the two ends come from the user input.
   - Masked-cell handling: exclude pixels of the mask color from the
     geometric median and mark fully masked cells as NaN.
4. [ ] New command-line driver `digitize_multi_image.py`: takes the image
   filename, approximate panels box, `n_rows`, `n_cols`, the per-row colorbar
   min/max value pairs, and output stem; writes per-panel outputs
   (`stem_r{r}_c{c}.txt/.npy`), a `grid.json` with all detected boxes and
   colorbar ends, and an optional debug overlay PNG.
5. [ ] Update `README.md` with the new script, its arguments, and an example
   based on the Arney figures.

Dependencies and order: steps 1 and 3 are independent and can be done in
parallel; step 4 needs 1 and 3; step 5 needs 4.

## Verification checklist

- `conda run -n digitize black --check` passes on all Python files.
- Regression: both `README.md` examples still reproduce `data_out.txt` and
  `data_out_cat.txt` exactly (single-image behavior must not change).
- Synthetic check: generate a fake panel-array figure with a known grid,
  colorbar, and mask fill, and verify the detected geometry and extracted
  values (including NaN placement).
- Real check on the Arney figures: inspect the debug overlay image
  (visually, by the user — agents cannot see images in this setup).

## Environment facts (verified)

- Conda environment for everything: `digitize`
  (`C:\Users\koeh_ti\AppData\Local\miniforge3\envs\digitize`). `conda` is not
  reachable inside the agent terminal sandbox; run conda commands unsandboxed.
- `black` 26.10.1 is installed in `digitize` (development-only).
- Working branch: `auto-multi-image` (first commit: `AGENTS.md`).
- `arney_fig26.png`, `arney_fig27.png`, and the Arney et al. 2014 PDF are
  untracked working-tree files (not in git, not ignored).

## Findings log

- (empty — to be filled during step 2; record verified pixel-level facts about
  the two figures here: sizes, background color, band structure, colorbar
  location, tick style)

## Resolved questions / gotchas

- `project_point_on_line_3d` contains a clamp `min(max(that, 0), -1)` that
  looks wrong but must not be changed as drive-by cleanup (see `AGENTS.md`).
- The categorical fallback in the original script hard-codes the output
  reshape to `(args.ny, args.nx)`; when refactoring, pass the shape instead.
- The current AI setup has no image vision: all figure inspection must be done
  numerically, and final visual checks are done by the user.
