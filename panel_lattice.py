"""
Define the panel lattice of a multi-panel figure from one reference box.

Everything about the figure reduces to four facts:

(a) a small reference box, taken inside the top-left panel, whose data area
    cross-correlates cleanly against the rest of the figure;
(b) the number of panel rows and columns;
(c) a coarse bounding box covering all panel data (the colorbar excluded);
(d) the position of every panel, measured as the shift of the reference box
    onto that panel and reported as absolute pixel coordinates relative to the
    top-left pixel of the image.

Because the panels form a rectangular grid, the reference box sits at the same
offset inside every panel. Measuring its absolute horizontal position along the
reference row and its vertical position along the reference column therefore
fixes all panels by their row and column. Each panel is larger than the
reference box; the height of its window is the equal-window size implied by the
coarse box, the row and column counts, and the white gaps between the panels,
while the width of a column window follows the white gutters of the panel
content: the boundary between two panel columns is a white gap of the whole
panel area, and the window edges sit at its centre, so every data pixel of a
column belongs to one window in every row even when the panels of a column
differ in how far sideways their data reaches. The
positions are measured by cross-correlating the signed (+1 data, -1 background
and masked pixels) field over the reference box, restricted to the expected
neighbourhood of each panel so that the strong repetition of the panels cannot
make a distant panel win a peak. The lattice is written to JSON, and the panel
windows, the crosses, and a translucent copy of the reference box are drawn
over the figure for a visual check.
"""

import argparse
import json

import numpy as np
from numpy.typing import NDArray
from scipy.signal import correlate

from digitize_gridded_image import load_image

# Manual lattice adjustments carried over from the coarse overlay the user
# accepted. The outer box is grown a little past the data box, and the equal
# window height separates the panel rows by a narrow white gap; the column
# gaps now only serve as the nominal guess for the anchor search. These
# constants reproduce the accepted geometry.
EXTRA_RIGHT_PX = 10
EXTRA_TOP_PX = 10
EXTRA_BOTTOM_PX = 10
GAP_X_PX = 10
GAP_Y_PX = 6

# The tolerance with which a pixel counts as white; the figures are clean
# digital renderings on a pure white background.
WHITE_TOLERANCE = 8


def runs(mask: NDArray[np.bool_]) -> list[tuple[int, int]]:
    """Inclusive index pairs of the contiguous True runs of a vector."""
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(np.int8), 0]))
    return list(zip(edges[0::2], edges[1::2] - 1))


def data_box(img: NDArray[np.int64]) -> tuple[int, int, int, int]:
    """
    The content box of the panel area, excluding the colorbar strip.

    Every non-white pixel outside the colorbars is data, so short slivers
    belong to the box as well; the colorbar strip is the right-most content
    column run of the figure.

    Parameters
    ----------
    img
        Image array with integer RGB channels.

    Returns
    -------
    box
        The (left, right, top, bottom) inclusive pixel bounds of the data.
    """
    nonwhite = np.any(np.abs(img - 255) > WHITE_TOLERANCE, axis=2)
    col_runs = [r for r in runs(nonwhite.any(axis=0)) if r[1] - r[0] + 1 >= 3]
    assert len(col_runs) >= 2, "expected panel content and a colorbar strip"
    left, right = col_runs[0][0], col_runs[-2][1]
    rows = np.flatnonzero(nonwhite[:, left : right + 1].any(axis=1))
    return int(left), int(right), int(rows[0]), int(rows[-1])


def column_windows(
    content: NDArray[np.float64],
    left: int,
    right: int,
    column_anchors: list[int],
    reference_width: int,
) -> tuple[list[int], list[int]]:
    """
    The window origin and width of every panel column, from white gutters.

    In the content projection of the panel area, the boundary between two
    neighbouring panel columns shows up as a white gap. The window of a
    column begins at the centre of the gap before it and ends at the centre
    of the gap after it, bounded by the outer box at the two ends, so that
    the shared gutter is split evenly between the neighbouring columns while
    every data pixel of each column stays covered. Inner white gaps of a
    panel (its fully masked cell columns) are excluded by restricting each
    boundary search to the region between the two measured column anchors.

    Parameters
    ----------
    content
        The 0/1 content mask of the whole figure, as from ``content_mask``.
    left
        The first column of the panel area (outer box left edge).
    right
        One past the last column of the panel area (outer box right edge).
    column_anchors
        The measured reference-box anchor position of every panel column, in
        image coordinates and in column order.
    reference_width
        The width of the reference box, which keeps an inner white gap of
        one panel from being mistaken for the gutter before the next.

    Returns
    -------
    origins
        The image ``x`` of the left edge of each panel column window.
    widths
        The width in pixels of each panel column window.
    """
    projection = (content[:, left:right] > 0).any(axis=0)
    white_runs = [(low + left, high + left) for low, high in runs(~projection)]
    centers: list[int] = []
    for anchor, next_anchor in zip(column_anchors, column_anchors[1:]):
        candidates = [
            gap
            for gap in white_runs
            if gap[0] >= anchor + reference_width and gap[1] < next_anchor
        ]
        assert candidates, f"no white gutter found between x {column_anchors}"
        gutter = max(candidates, key=lambda gap: gap[1] - gap[0])
        centers.append(int((gutter[0] + gutter[1] + 1) // 2))
    origins = [int(left)] + centers
    widths = [
        int(origins[index + 1] - origins[index])
        for index in range(len(column_anchors) - 1)
    ] + [int(right - origins[-1])]
    assert all(width > reference_width for width in widths), (
        f"column windows {widths} must all be wider than the reference box "
        f"({reference_width} px)"
    )
    return origins, widths


def content_mask(img: NDArray[np.int64]) -> NDArray[np.float64]:
    """
    The binary content mask of an image patch.

    The white background is False and every other colour is True, so the mask
    keeps the shape of the plotted area while ignoring its colours.

    Parameters
    ----------
    img
        Image patch with integer RGB channels.

    Returns
    -------
    mask
        Floating-point 0/1 mask, suitable for overlaying on the figure.
    """
    return np.any(np.abs(img - 255) > WHITE_TOLERANCE, axis=2).astype(np.float64)


def signed_mask(img: NDArray[np.int64]) -> NDArray[np.float64]:
    """
    The coloured/uncoloured field of an image patch for cross-correlation.

    Where the colour is that of the data (the RGB channels differ, so the pixel
    is neither white nor a shade of gray), the field is +1; everywhere else it
    is -1, marking both the white background and the masked black, gray, and
    enclosed-white pixels. Unlike the binary content mask, aligning two masked
    regions then counts as agreement, which is information about where the data
    is not.

    Parameters
    ----------
    img
        Image patch with integer RGB channels.

    Returns
    -------
    field
        Floating-point +1/-1 field, suitable for cross-correlation.
    """
    coloured = img.max(axis=2) - img.min(axis=2) > WHITE_TOLERANCE
    return np.where(coloured, 1.0, -1.0)


def best_position(
    profile: NDArray[np.float64], reference_size: int, nominal: float, search: int
) -> tuple[int, int, int]:
    """
    The best absolute coordinate near a nominal one, read from a profile.

    The search is confined around the nominal position so that the repetition
    of the figure cannot make a neighbouring panel win the peak, while still
    absorbing the drift of the true panels away from the equal-window guess.

    Parameters
    ----------
    profile
        The correlation over shift, indexed so that absolute coordinates are
        recovered from an index ``i`` as ``i - (reference_size - 1)``.
    reference_size
        The width or height of the reference box.
    nominal
        The nominal absolute coordinate to search around.
    search
        The half-width in pixels of the search window.

    Returns
    -------
    position
        The best absolute pixel coordinate, then the number of indices that
        attain the peak exactly and the number that reach 99% of it.
    """
    low = max(0, int(nominal - search) + reference_size - 1)
    high = int(nominal + search) + reference_size
    window = profile[low:high]
    peak = window.max()
    best = low + int(np.argmax(window))
    position = best - (reference_size - 1)
    exact = int(np.count_nonzero(window == peak))
    plateau = int(np.count_nonzero(window >= 0.99 * peak))
    return position, exact, plateau


def main() -> None:
    """Measure the panel lattice from the reference box and store it."""
    prs = argparse.ArgumentParser(description=__doc__ or "")
    prs.add_argument("filename", type=str, help="Filename of the figure image.")
    prs.add_argument("n_rows", type=int, help="Number of panel rows.")
    prs.add_argument("n_cols", type=int, help="Number of panel columns.")
    prs.add_argument(
        "anchor",
        type=int,
        nargs=4,
        metavar=("X", "Y", "WIDTH", "HEIGHT"),
        help="The reference box (a), inclusive top-left pixel and size, taken "
        "inside the top-left panel where the data cross-correlates cleanly.",
    )
    prs.add_argument(
        "--overlay",
        action="store_true",
        help="Draw the panel windows, the crosses, and the translucent "
        "reference box over the figure as a debug image.",
    )
    prs.add_argument(
        "--out-stem",
        type=str,
        default=None,
        help="Basename for the debug image and lattice file (defaults to the "
        "input filename).",
    )
    args = prs.parse_args()
    stem = args.out_stem or args.filename.rpartition(".")[0]
    ax, ay, aw, ah = args.anchor

    img = load_image(args.filename).astype(np.int64)
    field = signed_mask(img)
    content = content_mask(img)
    reference = field[ay : ay + ah, ax : ax + aw]

    left, right, top, bottom = data_box(img)
    # the outer box is the coarse data box grown by the accepted adjustments;
    # the growth at the top and left never runs off the image, whose pixels
    # do not exist - a figure whose data touches the top edge (its first row
    # of panel windows may start at y = 0) would otherwise give those windows
    # negative origins and empty image patches downstream; the growth at the
    # bottom and right may exceed the image, as consumers clip their patches
    outer_left = max(0, left)
    outer_right = right + 1 + EXTRA_RIGHT_PX
    outer_top = max(0, top - EXTRA_TOP_PX)
    outer_bottom = bottom + 1 + EXTRA_BOTTOM_PX
    # the equal-window sizes are the nominal guess of the anchor search; the
    # window height stays authoritative while the width is superseded by the
    # white gutters of the panel content below
    panel_w = (outer_right - outer_left - (args.n_cols - 1) * GAP_X_PX) / args.n_cols
    panel_h = (outer_bottom - outer_top - (args.n_rows - 1) * GAP_Y_PX) / args.n_rows
    print(
        f"outer box x {outer_left}..{outer_right - 1}, y {outer_top}..{outer_bottom - 1}; "
        f"panel height {panel_h:.2f} px, nominal width {panel_w:.2f} px, "
        f"gaps {GAP_X_PX}/{GAP_Y_PX}"
    )

    # the equal-window nominal position of the reference box corner in every
    # column and row, each refined to its measured absolute position
    nominal_x = [ax + col * (panel_w + GAP_X_PX) for col in range(args.n_cols)]
    nominal_y = [ay + row * (panel_h + GAP_Y_PX) for row in range(args.n_rows)]

    # One strip along the reference row gives every column anchor, and one
    # strip along the reference column gives every row anchor; each correlated
    # against the reference box a single time.
    strip_x = field[ay : ay + ah, :outer_right]
    profile_x = np.real(correlate(strip_x, reference, mode="full", method="fft"))[
        ah - 1, :
    ]
    columns: list[dict[str, int]] = []
    for col, nominal in enumerate(nominal_x):
        position, exact, plateau = best_position(
            profile_x, aw, nominal, int(panel_w // 2)
        )
        columns.append({"col": col, "x": position, "plateau_px": plateau})
        print(
            f"column {col}: x={position} (nominal {nominal:.0f}), "
            f"exact {exact}, 99% plateau {plateau} px"
        )

    strip_y = field[:outer_bottom, ax : ax + aw]
    profile_y = np.real(correlate(strip_y, reference, mode="full", method="fft"))[
        :, aw - 1
    ]
    rows_list: list[dict[str, int]] = []
    for row, nominal in enumerate(nominal_y):
        position, exact, plateau = best_position(
            profile_y, ah, nominal, int(panel_h // 2)
        )
        rows_list.append({"row": row, "y": position, "plateau_px": plateau})
        print(
            f"row {row}: y={position} (nominal {nominal:.0f}), "
            f"exact {exact}, 99% plateau {plateau} px"
        )

    # the reference box sits at a fixed offset inside each panel, so a panel
    # window of some row is its measured box shifted back by that offset; the
    # column windows themselves follow the white gutters of the panel area,
    # because panels of a column need not leave equal amounts of white at
    # their sides, and an anchor-derived width can cut off data or overlap
    column_x0_px, column_width_px = column_windows(
        content, outer_left, outer_right, [c["x"] for c in columns], aw
    )
    panel_y0 = [r["y"] - (ay - outer_top) for r in rows_list]
    print(f"column windows: origins {column_x0_px}, widths {column_width_px}")

    panels = [
        {"row": row, "col": col, "x": columns[col]["x"], "y": rows_list[row]["y"]}
        for row in range(args.n_rows)
        for col in range(args.n_cols)
    ]

    lattice = {
        "image": args.filename,
        "image_width_px": int(img.shape[1]),
        "image_height_px": int(img.shape[0]),
        "n_rows": args.n_rows,
        "n_cols": args.n_cols,
        "reference_box": {"x": ax, "y": ay, "width": aw, "height": ah},
        "outer_box": {
            "x0": int(outer_left),
            "x1": int(outer_right),
            "y0": int(outer_top),
            "y1": int(outer_bottom),
        },
        "panel_height_px": float(panel_h),
        "column_x0_px": column_x0_px,
        "column_width_px": column_width_px,
        "gap_px": {"between_columns": GAP_X_PX, "between_rows": GAP_Y_PX},
        "column_anchors_px": [c["x"] for c in columns],
        "row_anchors_px": [r["y"] for r in rows_list],
        "peak_plateau_px": {
            "columns": [c["plateau_px"] for c in columns],
            "rows": [r["plateau_px"] for r in rows_list],
        },
        "panels": panels,
    }
    json_name = f"{stem}_lattice.json"
    with open(json_name, "w", encoding="utf-8") as handle:
        json.dump(lattice, handle, indent=2)
        handle.write("\n")
    print(f"wrote {json_name}")

    if not args.overlay:
        return

    # delayed imports
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    reference_content = content[ay : ay + ah, ax : ax + aw]

    fig, view = plt.subplots(figsize=(img.shape[1] / 110, img.shape[0] / 110))
    view.imshow(
        img.astype(np.uint8),
        extent=(0, img.shape[1], img.shape[0], 0),
        interpolation="nearest",
    )

    windows = [
        [
            (x0, y0),
            (x0 + column_width_px[col], y0),
            (x0 + column_width_px[col], y0 + panel_h),
            (x0, y0 + panel_h),
            (x0, y0),
        ]
        for col, x0 in enumerate(column_x0_px)
        for y0 in panel_y0
    ]
    view.add_collection(LineCollection(windows, colors="black", linewidths=1.0))

    overlay = np.zeros((*img.shape[:2], 4))
    for panel in panels:
        x0, y0 = max(panel["x"], 0), max(panel["y"], 0)
        x1 = min(panel["x"] + aw, img.shape[1])
        y1 = min(panel["y"] + ah, img.shape[0])
        if x0 >= x1 or y0 >= y1:
            continue
        piece = (
            reference_content[
                y0 - panel["y"] : y1 - panel["y"], x0 - panel["x"] : x1 - panel["x"]
            ]
            > 0
        )
        overlay[y0:y1, x0:x1][piece] = (0.0, 0.0, 0.0, 0.2)
    view.imshow(
        overlay, extent=(0, img.shape[1], img.shape[0], 0), interpolation="nearest"
    )
    for panel in panels:
        view.plot(
            panel["x"],
            panel["y"],
            marker="x",
            color="black",
            markersize=12,
            markeredgewidth=2.0,
        )
    view.set_axis_off()
    filename = f"{stem}_lattice.png"
    fig.savefig(filename, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {filename}")


if __name__ == "__main__":
    main()
