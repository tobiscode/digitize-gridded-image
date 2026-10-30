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
  `n_rows` and `n_cols`, plus an anchor box (a small clean data window inside
  the top-left panel) from which the panel positions are cross-correlated
  (added 2026-10-30; see step 3 and the findings log).
- Output: full data extraction for every panel *and* the detected geometry
  (saved as JSON, plus an optional debug overlay image).
- Colorbar values: the current figures have *no tick marks*, so there is
  nothing to detect. Instead, the user provides the min/max values of each
  colorbar (one pair per row) separately; the script maps the two ends of the
  colorbar to those values. Which end holds min vs max must be explicit in the
  input (e.g. values given in top-to-bottom order) because the printed
  orientation is not readable from pixels; the debug output shows the color at
  each colorbar end so the user can sanity-check the orientation.
- Constraint set confirmed by the user on 2026-10-29 (supersedes earlier
  per-panel-phase and constant-gap reasoning):
  - ALL non-white pixels outside the colorbars are data (the figures are
    clean digital renderings, not scans), and NO data lies outside the
    lattice: the bounding box of the non-white pixels of the panel area IS
    the lattice extent, and every non-white pixel must be covered by a cell.
    A debug overlay showing non-white pixels outside the drawn lattice is a
    hard failure of the fit.
  - All cells have the same number of pixels in x, same number of pixels in y,
    same pixel width and same pixel height across the whole figure.
  - Coarse starting point: divide the non-white width and non-white height
    of the panel area by integer cell counts along each axis to get a first
    estimate of the individual cells.
- Panels are delimited only by white bands (no axis frames); under the
  zero-gap assumption those white bands are masked cell rows / columns.
- Masked cells (blank/fill-colored as well as black/grey cells) are reported
  as NaN.
- One colorbar at the far right of each row.

## Plan

Status markers: [ ] not started, [~] in progress, [x] done.

1. [x] Refactor `digitize_gridded_image.py` into importable functions:
   move the `__main__` pipeline into functions (image loading, per-cell
   dominant-color extraction, colorbar extraction, color-to-value mapping,
   saving, image recreation). Keep `geometric_median` and
   `project_point_on_line_3d` untouched (both were replaced by batched
   versions and deleted in the later vectorization, see the findings log),
   and keep the command-line interface
   (argument names and semantics) exactly as documented in `README.md`.
   (done, commit 5b0668242be2bb2a5600b834d6fc7ec681500ab9;
   outputs hash-identical to the reference files)
2. [x] Numerically verify the structure of `arney_fig26.png` /
   `arney_fig27.png` (dominant colors, white-gap row and column runs,
   colorbar location) to calibrate the detection parameters. See the
   Findings log for the measured values and detection parameters.
3. [x] New script `panel_lattice.py` (replaces the abandoned
   `grid_detection.py`, and the exploratory `coarse_lattice.py`,
   `panel_xcorr.py`, `content_lattice.py` - all deleted after the user
   approved the consolidated tool; they remain in git history). The panel
   lattice reduces to four facts: (a) a small reference anchor box inside the
   top-left panel where the data cross-correlates cleanly, (b) `n_rows` x
   `n_cols`, (c) the coarse data bounding box (colorbar excluded), and
   (d) the measured absolute position of the anchor box in each panel. The
   anchor box positions are measured by cross-correlating the SIGNED field
   (+1 coloured data, -1 background and masked) over the box only, along one
   row-strip and one column-strip, peak searched within half a panel of the
   equal-window nominal position; a rectangular grid then fixes all 36 panel
   positions, and panel windows get the equal-window height (accepted manual
   adjustments kept as constants: outer box +10 px right/top/bottom, gaps
   10 px x / 6 px y; the top/left growth is clamped to the image). The column
   windows follow the white gutters of the panel content projection (see
   `column_windows`): each boundary is the centre of the widest white gap
   between two neighbouring measured column anchors, so every data pixel of
   a column is covered in every row. Writes `<stem>_lattice.json` (anchors,
   column window origins and widths, panel height, gaps, per-panel boxes,
   peak plateau widths) and an optional `--overlay` debug PNG. Cell-size
   detection inside panels and value mapping are deliberately not in this
   script (they are steps 4 and 5).
4. [x] New scripts `cell_fft.py` + `cell_grid.py` fit the cell grid inside
   the panel windows (fig26 and fig27 measured 2026-10-30; both overlays
   accepted by the user). Panel window WIDTHS come from the white gutters of
   the column content (see step 3) - an equal-window width made columns 1/2
   of fig26 overlap and anchor-derived edges cut data off fig27; window
   heights are still the equal-window size. `cell_fft.py` stacks the
   |luminance-difference| amplitude spectra of all
   panel windows (band 3-9 px, parabolic peak) and stores the shared cell
   size in the lattice JSON (`cell_size_px`; renamed from `cell_pitch_px`,
   see the naming note below). `cell_grid.py` fits, per
   panel and axis, the phase of the integer seam lattice: nominal seams are
   `window origin + cell_size * k`, probed with sub-pixel
   shifts (0.02 px step) over one cell, each candidate scored by the MEAN
   absolute RGB jump across its seam lines (mean, not sum, so a phase that
   rounds an extra fringe line in gains nothing); reports the runner-up
   phase and per-line support as diagnostics, stores the integer seam lists
   per panel in the lattice JSON (`cell_grid`), and draws an `--overlay`
   debug PNG. Every panel keeps its own window-locked seam phase (JSON
   `cell_grid`; the seams of each figure also satisfy one global x lattice,
   see the findings log).
   Colorbar-end detection and value mapping happen in the driver of step 5.
5. [x] New command-line driver `digitize_multi_image.py` (done 2026-10-30;
   fig26/fig27 digitized and the recreations user-checked): takes the lattice
   JSON of steps 3-4, a bars JSON, and an output `.npy` file; digitizes every
   panel with the color logic of `digitize_gridded_image.py` and stacks the
   values into one array of shape (ny, nx, n_rows, n_cols), aligned by the
   seam lattice (whole-cell offsets of the panel seams from the earliest seam
   of its column or row band), with uncovered slots and masked cells NaN. The
   bars JSON names the color scale of each panel row,
   `{"bars": [{"label", "legend_box", "values"}, ...]}`, since the figures
   have no colorbar ticks; an optional `label` names its row in the
   recreation. Options: `--save-panels DIR` (per-panel `.txt` outputs),
   `--recreate-image FILE` (figure recreated panel by panel, each with its
   own small colorbar; drawn with GridSpec and constrained layout - which
   that layout only manages for `add_subplot`, not `add_axes`; saved at
   dpi*2), `--row-labels JSON` (else the bars labels are used),
   `--column-labels JSON`.
6. [x] Update `README.md` with the new script, its arguments, and an example
   (done 2026-10-30: the "Multi-panel figures" section documents the whole
   pipeline panel_lattice -> cell_fft -> cell_grid -> digitize_multi_image,
   including the bars JSON schema; it uses generic placeholder filenames -
   Arney data must not enter the repository).

Dependencies and order: steps 1 and 3 are independent and can be done in
parallel; step 4 needs 3; step 5 needs 1, 3 and 4; step 6 needs 5.

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
  (`C:\Users\koeh_ti\AppData\Local\miniforge3\envs\digitize`). Inside the
  agent terminal sandbox the environment's `python.exe` runs when PATH is
  prepended with `envs\digitize`, its `Scripts`, and its `Library\bin`
  (miniforge3 is sandbox-readable via the user's
  `chat.agent.sandbox.fileSystem.windows.allowRead` setting); `conda.exe`
  itself cannot start in the sandbox (system-DLL relocation error), and
  `conda run` also withholds the child's output until it exits, which makes
  logs look dead - prefer `python.exe -u` directly.
- `black` 26.10.1 is installed in `digitize` (development-only). It HANGS in
  the sandbox (its cache/temp writes are blocked); run it unsandboxed (via
  `conda run -n digitize black`).
- Arney data must NOT enter the repository: all `arney_*` files (figures,
  lattices, bars, outputs) stay untracked, README examples use generic
  placeholders, and the user performs all commits - agents stage nothing.
- Working branch: `auto-multi-image` (first commit: `AGENTS.md`).
- `arney_fig26.png`, `arney_fig27.png`, and the Arney et al. 2014 PDF are
  untracked working-tree files (not in git, not ignored).

## Findings log

(verified numerically on 2026-10-29 with temporary probe scripts, since there
is no agent image vision; user confirmed there is NO text in either figure)

- `arney_fig26.png` is 1351x2352, `arney_fig27.png` is 1368x1691. Background
  is pure white (255,255,255); a channel tolerance of 8 reproduces the
  non-background mask cleanly. White also doubles as the mask-fill color.
- Lattice via projection runs of the non-background mask:
  - fig26: exactly 9 row bands (starts 1, 265, 529, 793, 1057, 1322, 1586,
    1851, 2115; pitch ~264 px, gaps >= 4 px). Panel column groups at
    x ~ 27-220, 302-512, 572-825, 891-1145; colorbar strip at x ~ 1293-1350.
  - fig27: exactly 9 row bands (0..1690, pitch ~188 px); 7 panel column
    groups at x ~ 0-98, 166-290, 349-482, 557-674, 717-866, 925-1058,
    1117-1266; colorbar strip at x ~ 1320-1367. A morphological opening
    (binary erosion 2 px then dilation 2 px) is needed before projections:
    without it, panels/bar segments of adjacent rows merge vertically.
- The column "fragments" separated by >= 15 px white gaps within the panel
  area are fully masked cell columns of the panels (NOT text), e.g. fig26
  col group starts at x=27 then a gap before x=93. Panel rectangles are
  therefore the column-group spans, interior white strips included.
- Cell pitch (autocorrelation of pooled cell-seam profiles; profile per
  internal boundary column/row using mean |dRGB| over pairs where both sides
  are non-white):
  - fig26: fundamental ~5.5 px in both directions (harmonic chain 5-6, 11,
    16-17, 22, 27-28; fundamental confirmed by per-panel seam spacings).
  - fig27: exactly 4 px (chain 4, 8, 12, 16, 20 across all column groups AND
    the row direction).
  - Seams are 1-px antialiased lines; masked regions intrude INTO cells
    (mask covers one side of a cell), shifting naive dominant colors -
    near-white pixels must be excluded per cell; fully masked cells -> NaN.
- Panel pixel widths DIFFER between column groups: fig26 194/211/254/255 px
  with common pitch 5.54; fig27 ~96/133/134/118/150/134/150 px at pitch 4.
  Under the zero-gap model those groups are separated by columns of fully
  masked cells, and every panel still holds the same `nx`.
- Colorbars: continuous vertical gradient, no ticks, 58 px (fig26) / 44 px
  (fig27) wide; ~387 unique colors per column in fig26. Bar vertical extent
  differs from the panel extent of the same row band (fig26 row 0: bar rows
  1-238, panel rows 33-259), so the bar ends must be measured inside the bar
  columns, not copied from the band. fig27's colorline is scan-noisy
  (median adjacent-pixel |dRGB| = 3) -> take a column median across the bar's
  interior columns plus a small median filter along the bar.
- 2 user-supplied end values per colorbar would make only 2 markers; the
  existing `build_index_to_value` CubicSpline needs >= 3 points and would
  silently fall back to categorical - fix by passing >= 4 collinear pseudo-
  markers (linear ramp between the user's min/max), which reproduces a
  exactly linear map through the existing spline code path.
- Panel-area non-white box (the lattice extent per the user's constraints):
  fig26 x 0..1145, y 33..2351; fig27 x 0..1263, y 0..1690 (colorbars
  excluded). The short black bar at x=0..5 in rows 254-259 and 515-520 of
  fig26 is DATA belonging to the left-most cell columns (user: no scan
  artifacts), so the lattice's left edge is x=0, not the first visible seam
  near x=27.
- The seam-F-score pitch fits of the superseded `grid_detection.py` (shared
  pitch with per-panel phases: fig26 5.55/5.00 px, fig27 4.10/3.90 px) were
  ABANDONED: the user's debug overlay showed data pixels outside the drawn
  lattice. The later zero-gap box-tiling model was likewise abandoned because
  geometry alone cannot fix the panel grid.
- FINAL anchor-box method for fig26 (verified 2026-10-30, user approved the
  overlay): anchor box x=93, y=33, 104x221 px inside panel (0,0) (the left
  edge of its semicircle is the one feature consistent across panels); the
  ±1 signed-field correlation over that box measures every panel anchor
  within ~1 px of the earlier manually locked grid. Measured anchors:
  column x = 93, 380, 655, 974 (equal-window nominal 93, 384, 676, 968);
  row y = 33, 295, 556, 818, 1079, 1341, 1603, 1865, 2126. All column peaks
  and row peaks are 1-px-wide exact maxima (99% plateau <= 8 px).
- Lessons for the column x-positions of binary/overlap correlation (why the
  plain route failed): each panel shows a semicircle with a different amount
  cut away (content widths ~221-270 px), so a max-overlap correlation
  SATURATES - the wider mask covers the narrower over a range of shifts,
  giving 13-28 px flat plateaus and meaningless peaks (measured, not
  assumed). Fix that worked: correlate only a small clean data box (the
  anchor) instead of the whole panel, AND use the signed field so agreeing
  masked pixels count; this collapsed the plateaus to <= 11 px and moved
  column 2 off its wrong zero. Far-left scattered black/gray pixels at
  x=0..5, rows 254-259/515-520 are NOISE for alignment purposes (user,
  2026-10-30), even though they are data for the box-extent purpose.
- Cell pitch (stacked FFT of the panel difference profiles, re-measured on
  the gutter windows 2026-10-30; the last digit follows the window crop): fig26
  w_x = 5.5104 px, w_y = 5.0032 px; fig27 w_x = 4.0917, w_y = 3.9261 px. The
  per-panel spectral peaks agree within ~0.01 px per axis, so one shared pitch
  per figure is justified.
- Seam-phase fit of fig26 (verified 2026-10-30, user accepted the overlay):
  the fitted window-relative shift is nearly constant per panel, y phase
  alternating 0.40/0.48 px by row; the x phase reads from each gutter
  window's origin. Keep per-panel seam lists. Open anomaly, deliberately not
  tie-broken yet: panels (4,2)-(4,3) choose y shift 1.40 (51 seams) against
  the neighbor phase 0.48 even under mean scoring, with normal line support -
  row 4 may genuinely sit ~1 px lower in its windows. (4,1) chose 0.48 once
  the column windows stopped overlapping into column 1.
- Panel window WIDTHS must be per column (verified 2026-10-30, user saw the
  overlap on the equal-window overlay): the measured column anchor steps are
  287/275/319 px, so an equal window of 281.5 px overlapped columns 1/2 by
  6.5 px and left a 37.5 px air gap at column 2/3. Windows now run anchor to
  next anchor minus GAP_X_PX (widths 277/265/309/275 px, last column to the
  outer box right edge); a coverage probe found 0 non-white data pixels
  outside them - the equal-window width, not the figure, was breaking the
  global x lattice. Row anchor steps 261-262 px do NOT violate the
  equal-window height (254.56 + 6 gap), so per-row heights are not needed.
- Panel window LEFT edges come from the white gutters of the content
  projection, not from the anchor (verified 2026-10-30, user asked for it
  after the fig27 coverage probe): anchor-derived origins cut data off,
  because panels in one column can differ in how far their data reaches
  sideways (fig27 lost 59 245 px, up to 48 px = ~12 cell columns of the
  narrow last column). A gutter is the widest white gap between two
  neighbouring measured column anchors, and the window boundary is its
  centre; inner white gaps of a panel (fully masked cell columns) are
  excluded by that anchor restriction. fig26 gutter centres land exactly on
  the true content gaps (origins 0/261/542/858, widths 261/281/316/298) and
  it keeps 0 uncovered pixels; fig27 origins 0/128/320/504/696/896/1088,
  widths 128/192/184/192/200/192/186. The x seam phase is now nominal from
  the window origin (0.0), not from the anchor offset, since gutters hold no
  fixed offset inside a panel.
- One global x seam lattice DOES explain both figures, which the per-panel
  phase metric hides (least-squares fit of all vertical seams in image
  coordinates, 2026-10-30): fig26 1872 seams fit pitch 5.5025 with residuals
  rms 0.34 / max 0.96 px; fig27 2781 seams fit pitch 4.0907 with residuals rms
  0.34 / max 0.83 px. The gutter origins are not round lattice multiples, so
  the circular per-panel phase agreement reads coarsely (fig26 1.14 px, fig27
  0.63 px x) while the seams themselves stay put. The cell grid keeps the
  per-panel integer seam lists either way.
- fig27 lattice (measured 2026-10-30, anchor box x=4,y=20,85x156 in panel
  (0,0) - a box inside a middle column is NOT allowed, panel_lattice nominalx
  assumes the reference box starts in column 0 and the search window ran off
  the image): column anchors x = 4,197,384,575,768,965,1172, row anchors
  y = 20,208,...,1529 (steps 188 px, all 1-px exact peaks). fig27 data touches
  the image top, so the outer box top growth clamps to y=0 (see the clamp in
  main). After the gutter fix 473 px of fig27 data still fall outside the row
  window bands, all 1-3 px beyond a band edge (panels bleeding into the 6 px
  row gaps) - a row-direction matter, separate from the column gutters.
- Convention of the seam positions (measured on raw pixels, not guessed):
  a seam at position s scores the pixel PAIR (s-1, s); the visible 1-px
  antialiased blend pixel is row/column s-1, so debug lines must be drawn
  through coordinate s-0.5 (imshow pixel convention, where pixel r spans
  [r, r+1]). The digitizer must sample cell interiors between seam blend
  pixels accordingly. Sanity check of the fit: mean jump of the chosen pair
  2.65 (vertical seams) / 3.00 (horizontal) versus 0.5-1.0 one px away in
  either direction.
- Naming of the cell size (user choice, 2026-10-30): the shared cell size is
  `cell_width`/`cell_height` in the code and the JSON key `cell_size_px`
  (top level of the lattice file, and inside the `cell_grid` block); the old
  `pitch` wording (`cell_pitch_px`, `pitch_px`) is gone from all scripts and
  from both lattice JSONs, which were regenerated and verified numerically
  identical apart from the key rename.
- Batched dominant-color extraction (verified 2026-10-30, both README
  examples byte-identical after each step): `dominant_cell_colors` runs ONE
  batched `geometric_medians` call per panel (per-pixel owner arrays over
  the whole window, seam pixels belong to no cell, cells without a coloured
  pixel are masked), and `extract_dominant_colors` is likewise loop-free. The
  copied scalar `geometric_median` was deleted with its `cdist`/`euclidean`
  imports, and its attribution comment now sits above `geometric_medians`
  (keep it, see `AGENTS.md`). `project_points_on_line_3d` was vectorized in
  the same way (2.27 s -> 60 ms per panel) and its scalar twin deleted.
- Dominant colors are ROUNDED, not truncated: both paths end in
  `.clip(0, 255).round().astype(int)`. Bare `astype(int)` truncation turned
  Weiszfeld float noise (~1e-6) into off-by-one colors with a systematic
  DOWNWARD bias on flat scan cells. Batch versus old scalar paths still
  differ by +-1 on a few real-scan cells (different summation order, and
  geometric medians of degenerate cell sets are not unique); accepted, the
  masked sets always stay identical.
- Weiszfeld iteration cap (added 2026-10-30): a batch stops only when every
  cell converged, and a few flat-colour-plus-scan-noise cells creep
  asymptotically and stall their whole batch (~10 s/panel uncapped). In
  `geometric_medians`, `max_iter=1000` caps the passes; cells still moving
  then keep their latest estimate (movement is far below a color unit).
  Every call prints "N passes, X of Y cells converged by themselves[, Z
  stopped at the cap of 1000]". Measured: the README examples need 163
  (continuous) and 27 (categorical) passes, so all four reference files
  stay byte-identical; capped fig26/fig27 runs need ~40 s each, with 1.54 %
  of 32 435 cells (worst panel 3.6 %) and 0.77 % of 58 746 cells (worst
  1.6 %) stopping at the cap.

## Resolved questions / gotchas

- `project_points_on_line_3d` contains a clamp `min(max(that, 0), -1)` that
  looks inverted but is deliberate and passes NaN through unchanged; it must
  not be "fixed" as drive-by cleanup (see `AGENTS.md`). Its scalar-twin
  wrapper no longer exists.
- The categorical fallback in the original script hard-codes the output
  reshape to `(args.ny, args.nx)`; when refactoring, pass the shape instead.
- The current AI setup has no image vision: all figure inspection must be done
  numerically (agents draw overlay PNGs and print tables), and final visual
  checks are done by the user. The user's workflow rule: one small step per
  request, ask follow-up questions, no unrequested extra analysis.
- Run regression checks with the EXACT commands from `README.md`; once,
  reconstructing the categorical example's arguments from memory produced a
  false "regression" (the real categorical command is the long one there,
  with box [117, 379, 313, 574] and labels A-G).
