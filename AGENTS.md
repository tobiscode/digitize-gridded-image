# AGENTS.md — Digitize Gridded Images

Instructions for AI coding agents working in this repository.

## Purpose of the project

This repository contains a tool for digitizing gridded (raster) scientific plots:
given an image whose colored patches were generated from a 2D data field and the
matching colorbar, it recovers the underlying numerical values from the colors.
The main script is `digitize_gridded_image.py`. Audiences are scientists (e.g.
earth and climate scientists) who use the tool from the command line, not software
engineers.

## Environment

- Use the conda environment **`digitize`** for *everything* (running the script,
  tests, formatting, inspections). Never assume the base environment.
- Invoke it non-interactively, e.g.:

    ```powershell
    conda run -n digitize python digitize_gridded_image.py ...
    conda run -n digitize black <files>
    ```

- Dependencies are pinned in `requirements.txt` (numpy, pillow, scipy,
  matplotlib). Do not add new runtime dependencies without explicit confirmation;
  if a development-only tool is needed (e.g. `black`), install it into `digitize`
  and mention it rather than silently adding it to `requirements.txt`.
- `black` is installed in `digitize` as a development-only tool (not part of
  `requirements.txt`).

## Branch

- Active development happens on the branch **`auto-multi-image`** (aim: automatic
  handling of multiple images; `AGENTS.md` is its first commit). Do the work
  there rather than on `main`.

## Code style

- Format all Python code with **`black`** (default settings, 88 columns). Run
  `conda run -n digitize black <file>` after every edit.
- **Always use type hints** on every function signature (parameters and return
  types), using modern Python 3.11+ syntax: builtin generics (`list[int]`,
  `tuple[int, np.ndarray]`) and `X | None` instead of `typing.Optional`. Use
  `numpy.typing.NDArray` or precise `np.ndarray` annotations where helpful.
- Document functions, classes, and modules with **numpydoc-style docstrings**:
  one-line summary, then sections such as `Parameters`, `Returns`, `Raises`,
  `Notes`, `Examples` as appropriate. Document types in docstrings as well, in
  numpydoc style (e.g. `ndarray`), consistent with the hints.
- Keep terminology scientific rather than software-engineering oriented:
  e.g. "data field", "grid cell", "colorbar", "pixel coordinates", "RGB space",
  "dominant color" — avoid jargon like "payload", "handler", "util", "DTO".
  Variable names should follow the same convention (`colorline`, `markers_data`,
  `borders_x`), consistent with the existing code.
- Comments and docstrings are written in plain English, complete sentences,
  explaining *why* an algorithmic step is done (existing style).

## Structure

- Refactoring into importable functions and modules is allowed and encouraged
  where it improves clarity, but the command-line workflow described in
  `README.md` must keep working with the same argument names and semantics.
- If the CLI or any documented behavior changes, update `README.md` in the same
  change, including its example invocations.

## Validation (required after every change)

1. Formatting:

    ```powershell
    conda run -n digitize black --check digitize_gridded_image.py
    ```

2. Functional check — run the README example end to end and compare against the
   committed reference outputs (`data_out.txt` for the continuous example,
   `data_out_cat.txt` for the categorical one):

    ```powershell
    conda run -n digitize python digitize_gridded_image.py image_in.png "[498, 379, 693, 574]" "[715, 377, 715, 577]" "[377, 427, 477, 527, 577]" "[1, 0.5, 0, -0.5, -1]" 7 7 data_out.txt --recreate-image=image_out.png
    ```

   There is no formal test suite; this README example is the de-facto test. If
   `image_in.png` is missing from the workspace, obtain or recreate it as
   described in `README.md` before validating, or tell the user it is needed.

## Notes on the existing implementation

- `geometric_median` is copied code (StackOverflow answer, CC-BY-SA 3.0 / zlib);
  keep its license comment intact.
- `project_point_on_line_3d` contains a deliberate-looking clamp
  `min(max(that, 0), -1)`; verify against tests/examples before "fixing" it, and
  do not change it as drive-by cleanup.
- Colors are handled as integer RGB triples in the range 0–255; the colorbar is
  assumed to be a perfectly vertical or horizontal line in pixel space.
