"""
Digitize every panel of a multi-panel figure into one four-dimensional array.

The panel lattice (``panel_lattice.py``) and the cell grids inside the panel
windows (``cell_fft.py`` and ``cell_grid.py``) are taken as read from the
lattice JSON. This script then applies the logic of
``digitize_gridded_image.py`` to every panel: the dominant colour of every cell
between neighbouring seams is extracted, located on the colorbar of the panel
row, and converted to a data value using the two end values the user supplies
for that colorbar (the figures have no tick marks). The extraction of the
original script is reused unchanged; only the grid of cells comes from the
fitted seams instead of an evenly divided box.

Because the panels of a column share their vertical seams and the panels of a
row share one colorbar, the values of all panels stack into one array of shape
(ny, nx, n_rows, n_columns): nx is the largest cell count reached by any
column window and ny that of any row band, aligned by the seam lattice (whole
cell offsets of the panel seams from the earliest seam of its column or row
band), with uncovered slots and masked cells set to NaN. A cell is masked when
it holds no coloured pixel at all, which covers the white mask fill, the
enclosed white, and black or gray mask shades alike.

With ``--recreate-image``, the stacked array is drawn back as a picture, in the
spirit of ``digitize_gridded_image.recreate_image``. The recreation does not
know that the panels of a row share their colorbar - in general every panel
has its own - so every rendered panel gets its own small colorbar next to it,
even though the panels of a row are colored with the same scale.
"""

import argparse
import json
import os
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from digitize_gridded_image import (
    RGBImage,
    build_index_to_value,
    colors_to_colorline_indices,
    geometric_medians,
    load_image,
    map_indices_to_values,
    save_values,
    unique_colors,
)
from panel_lattice import WHITE_TOLERANCE

if TYPE_CHECKING:
    from matplotlib.figure import Figure

# the number of evenly spaced pseudo-markers placed along a colorbar between
# its two user-supplied end values; the original script needs at least three
# markers to build a continuous (cubic, here exactly linear) value map, and
# two values alone would silently fall back to categorical labels
N_PSEUDO_MARKERS = 5


def dominant_cell_colors(
    img: RGBImage, borders_x: NDArray[np.int64], borders_y: NDArray[np.int64]
) -> tuple[NDArray[np.int64], NDArray[np.bool_]]:
    """
    The dominant colour of every cell, and which cells are masked out.

    The cell between two consecutive seams takes the pixels strictly between
    them, excluding the antialiased seam pixel of the lower seam (a seam at
    position ``s`` scores the pixel pair (s - 1, s), whose visible blend pixel
    is row/column s - 1). A pixel counts as coloured when its RGB channels
    differ by more than the white tolerance, the data/background distinction
    of ``panel_lattice.signed_mask``, so white mask fill, background, enclosed
    white, and black or gray mask shades are all uncoloured. A cell with no
    coloured pixel is masked: its colour estimate is left at zero and it cannot
    be located on a colorbar meaningfully. Otherwise its colour is the
    geometric median of the cell's coloured pixels as in
    ``digitize_gridded_image.extract_dominant_colors``.

    Parameters
    ----------
    img
        Image array with integer RGB channels.
    borders_x
        The absolute seam positions delimiting the cells horizontally.
    borders_y
        The absolute seam positions delimiting the cells vertically.

    Returns
    -------
    colors
        Array of shape (ny, nx, 3) with the dominant integer RGB colour of
        every cell, zero in masked cells.
    masked
        Boolean array of shape (ny, nx), True where the cell was masked.
    """
    ny = len(borders_y) - 1
    nx = len(borders_x) - 1
    # the cell of every row and column of the window, -1 for the antialiased
    # seam pixel, which belongs to no cell and is skipped, mirroring the
    # `borders[i] : borders[i + 1] - 1` slice of the per-cell extraction
    rows = np.arange(borders_y[0], min(borders_y[-1] - 1, img.shape[0]))
    columns = np.arange(borders_x[0], min(borders_x[-1] - 1, img.shape[1]))
    if rows.size == 0 or columns.size == 0:
        # a window reaching past the image leaves no pixel for any cell
        return (
            np.zeros((ny, nx, 3), dtype=np.int64),
            np.ones((ny, nx), dtype=bool),
        )
    owner_y = np.searchsorted(borders_y, rows, side="right") - 1
    owner_x = np.searchsorted(borders_x, columns, side="right") - 1
    owner_y = np.where(rows == borders_y[owner_y + 1] - 1, -1, owner_y)
    owner_x = np.where(columns == borders_x[owner_x + 1] - 1, -1, owner_x)
    block = img[rows[0] : rows[-1] + 1, columns[0] : columns[-1] + 1, :]
    inside = (owner_y[:, None] >= 0) & (owner_x[None, :] >= 0)
    pixels = block[inside]
    owner = (
        np.broadcast_to(owner_y[:, None], inside.shape)[inside] * nx
        + np.broadcast_to(owner_x[None, :], inside.shape)[inside]
    )
    # a pixel counts as coloured when its RGB channels differ by more than the
    # white tolerance, the data/background distinction of
    # ``panel_lattice.signed_mask``, so white mask fill, background, enclosed
    # white, and black or gray mask shades are all uncoloured; only coloured
    # pixels take part in the median of their cell
    coloured = np.ptp(pixels, axis=1) > WHITE_TOLERANCE
    pixels, owner = pixels[coloured], owner[coloured]
    # the dominant colour of every cell at once, and which cells hold no
    # coloured pixel at all and are therefore masked
    medians = geometric_medians(pixels, owner, ny * nx)
    # round rather than truncate the medians: the iterative median carries
    # float noise far below a color unit which bare truncation would turn
    # into off-by-one colors whenever the true value is an integer
    colors = np.clip(medians, 0, 255).round().astype(np.int64).reshape(ny, nx, 3)
    counts = np.bincount(owner, minlength=ny * nx)
    masked = (counts == 0).reshape(ny, nx)
    return colors, masked


def colorline_of_bar(img: RGBImage, legend_box: list[int]) -> NDArray[np.int64]:
    """
    The colors along a colorbar, robust against noise and edge bleed.

    The original ``extract_colorline`` reads the single pixel line of the
    legend box, which suffers both from the scan noise of a printed figure and
    from the contaminated first and last pixel of the bar. Instead, every row
    of the bar is taken as the median across the bar strip around the supplied
    ``x``, the sequence is median-filtered along its length, and its first and
    last pixel are dropped. The colorbar is assumed to be vertical, as in
    these figures.

    Parameters
    ----------
    img
        Image array with integer RGB channels.
    legend_box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        colorbar as in ``digitize_gridded_image.py``, with left = right.

    Returns
    -------
    colorline
        Array of shape (length, 3) with the integer RGB color of every pixel
        along the colorbar from top to bottom, inset by one pixel at each end.
    """
    box = np.asarray(legend_box)
    assert box[0] == box[2], "this script expects vertical colorbars"
    assert box[3] > box[1] + 2, "the colorbar is too short to inset its ends"
    column = int(round(box[0]))
    half_width = 4
    strip = img[
        box[1] : box[3] + 1,
        max(0, column - half_width) : column + half_width + 1,
        :,
    ].astype(np.float64)
    colorline = np.median(strip, axis=1)
    filtered = colorline.copy()
    for index in range(1, len(colorline) - 1):
        filtered[index] = np.median(
            np.stack([colorline[index - 1], colorline[index], colorline[index + 1]]),
            axis=0,
        )
    return filtered[1:-1].round().astype(np.int64)


def pseudo_markers(
    n_colorline: int, values: list[float]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Collinear pseudo-markers turning the two end values into a linear map.

    Parameters
    ----------
    n_colorline
        Number of pixels of the extracted colorbar.
    values
        The data values [top, bottom] at the two ends of the colorbar.

    Returns
    -------
    positions
        Ascending positions evenly along the colorbar.
    data
        The linear ramp of data values between the supplied ends.
    """
    positions = np.linspace(0.0, n_colorline - 1, N_PSEUDO_MARKERS)
    data = np.linspace(float(values[0]), float(values[1]), N_PSEUDO_MARKERS)
    return positions, data


def extract_panel(
    img: RGBImage,
    bar: dict[str, object],
    borders_x: NDArray[np.int64],
    borders_y: NDArray[np.int64],
) -> tuple[NDArray[np.float64], int]:
    """
    The data values of one panel, and the number of masked cells.

    Parameters
    ----------
    img
        Image array with integer RGB channels.
    bar
        The colorbar entry of the panel row: the legend box of its colorbar
        and the data values [top, bottom] at its two ends.
    borders_x
        Absolute seam positions delimiting the cells horizontally.
    borders_y
        Absolute seam positions delimiting the cells vertically.

    Returns
    -------
    values
        Array of shape (ny, nx) with the data value of every cell and NaN in
        the masked cells.
    n_masked
        The number of masked cells.
    """
    colors, masked = dominant_cell_colors(img, borders_x, borders_y)
    colorline = colorline_of_bar(img, bar["legend_box"])  # type: ignore[arg-type]
    positions, data = pseudo_markers(len(colorline), bar["values"])  # type: ignore[arg-type]
    index_to_value, categorical = build_index_to_value(positions, data)
    assert not categorical, f"colorbar {bar['legend_box']} mapped to categorical"
    _, indices_general = colors_to_colorline_indices(
        colors, colorline, unique_colors(colorline)
    )
    values = map_indices_to_values(
        index_to_value, False, indices_general, indices_general.shape
    ).astype(np.float64)
    values[masked] = np.nan
    return values, int(masked.sum())


def window_borders(
    panel_grid: dict[str, object],
    axis: int,
    window_start: int,
    window_end: int,
    cell_size: float,
) -> NDArray[np.int64]:
    """
    The absolute seam positions of one panel axis, extended over the window.

    Cells live between consecutive seams; the fragments of the outermost cells
    before the first and after the last fitted seam are recovered by
    continuing the seam lattice by whole cells, so that no pixel of a panel
    is left unassigned. Fragments reaching into the white border are masked
    out as ordinary uncoloured cells.

    Parameters
    ----------
    panel_grid
        One entry of the lattice JSON ``cell_grid`` of the panel.
    axis
        0 for the rows (y) and 1 for the columns (x) of the cells.
    window_start
        The absolute image coordinate of the first window pixel.
    window_end
        The absolute image coordinate one past the last window pixel.
    cell_size
        The shared cell size along this axis in pixels, its width for the
        columns and its height for the rows.

    Returns
    -------
    borders
        The absolute seam positions delimiting the cells of the panel.
    """
    key = "seams_x" if axis == 1 else "seams_y"
    origin = int(panel_grid["x0"]) if axis == 1 else int(panel_grid["y0"])
    seams = [origin + seam for seam in panel_grid[key]]  # type: ignore[operator]
    low = seams[0]
    while low - cell_size >= window_start:
        low = int(round(low - cell_size))
        seams.insert(0, low)
    high = seams[-1]
    while high + cell_size <= window_end:
        high = int(round(high + cell_size))
        seams.append(high)
    return np.array(seams, dtype=np.int64)


def panel_placement(
    borders_x: NDArray[np.int64],
    borders_y: NDArray[np.int64],
    origin_x: int,
    origin_y: int,
    cell_width: float,
    cell_height: float,
    nx_big: int,
    ny_big: int,
) -> tuple[int, int, int, int]:
    """
    Where a panel's cells land in the stacked array, and how many fit.

    A panel is placed by the whole-cell offset of its first seam from the
    earliest seam of its column or row band, which aligns the panels at their
    seam lattice rather than at their cross-correlation anchors; panels that
    fit one cell more (their lattice is displaced by one cell) take one extra
    slot, as agreed.

    Parameters
    ----------
    borders_x
        Absolute seam positions of the panel horizontally.
    borders_y
        Absolute seam positions of the panel vertically.
    origin_x
        Absolute position of the earliest seam of the panel's column.
    origin_y
        Absolute position of the earliest seam of the panel's row band.
    cell_width
        Shared cell width in pixels.
    cell_height
        Shared cell height in pixels.
    nx_big
        Number of stacked columns available.
    ny_big
        Number of stacked rows available.

    Returns
    -------
    placement
        The (column, row) slot of the panel's first cell in the stack, then
        the number of its cells that fit horizontally and vertically, clamped
        to the stack.
    """
    offset_x = int(round((borders_x[0] - origin_x) / cell_width))
    offset_y = int(round((borders_y[0] - origin_y) / cell_height))
    nx = min(len(borders_x) - 1, nx_big) - offset_x
    ny = min(len(borders_y) - 1, ny_big) - offset_y
    return offset_x, offset_y, max(0, nx), max(0, ny)


def check_bar(bar: dict[str, object]) -> None:
    """
    Validate one colorbar entry of the bars file.

    Parameters
    ----------
    bar
        The entry with the legend box of the colorbar and the data values at
        its two ends.

    Raises
    ------
    AssertionError
        If the legend box is not a vertical line or the values are unusable.
    """
    box = np.asarray(bar["legend_box"])
    assert box.size == 4, "a legend box needs [left, top, right, bottom]"
    assert box[0] == box[2], "the legend box has to be a vertical line"
    assert box[3] > box[1], "the legend box is empty"
    values = np.asarray(bar["values"], dtype=float)
    assert values.size == 2, "two values are needed, one per colorbar end"
    assert np.all(np.isfinite(values)), "the colorbar values have to be numbers"
    assert not np.isclose(values[0], values[1]), "the two colorbar values are equal"


def recreate_figure(
    stack: NDArray[np.float64],
    boxes: list[tuple[int, int, int, int, int, int]],
    bars: list[dict[str, object]],
    colorlines: list[NDArray[np.int64]],
    cell_width: float,
    cell_height: float,
    row_labels: list[str] | None,
    column_labels: list[str] | None,
    title: str,
    filename: str,
    dpi: float,
    save_scale: float = 2.0,
) -> "Figure":
    """
    Draw the stacked array back as a figure, one small colorbar per panel.

    Every panel is rendered into an axes placed and sized by its slot in the
    stacked array, so that the picture reproduces the geometry of the original
    figure; the cells are drawn as colored squares, and the missing cells
    (masked and uncovered alike) stay blank. Since the recreation does not know
    that the panels of a row share their colorbar, every panel is given its own
    colorbar, like in the original figures - the panels of a row nevertheless
    are colored by the same sequence of colors and the same two end values,
    because those are what their shared colorbar provides.

    Parameters
    ----------
    stack
        The values of shape (ny, nx, n_rows, n_cols) as saved by this script.
    boxes
        One ``(row, col, offset_y, offset_x, ny, nx)`` tuple per rendered
        panel, its slot in the stacked array, as placed while digitizing.
    bars
        The colorbar entries of the panel rows, as read from the bars file,
        providing the two end values of the color scale of each row.
    colorlines
        The colors of the colorbar of every panel row, as extracted by
        `colorline_of_bar`, from its top end to its bottom end.
    cell_width
        Shared cell width in pixels, to scale the panels.
    cell_height
        Shared cell height in pixels, to scale the panels.
    row_labels
        The label of every panel row, or None to leave the rows unnamed.
    column_labels
        The label of every panel column, or None to leave the columns unnamed.
    title
        A title to draw above the figure, possibly empty.
    filename
        Where to save the recreated figure.
    dpi
        Resolution of the saved figure, in pixels per inch.
    save_scale
        Factor by which the figure is saved above its layout resolution, so
        that the small print of the many colorbars stays readable.

    Returns
    -------
    figure
        The drawn figure, already saved, ready for inspection.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    n_rows, n_cols = stack.shape[2], stack.shape[3]

    left_in = 1.8 if row_labels else 0.4
    top_in = 0.9 if title else 0.3
    bottom_in = 0.5 if column_labels else 0.3
    # the slot of a column or row band is as large as the cells its panels
    # occupy, measured from the far edge of the last cell placed in it
    slot_w = [
        max((box[3] + box[5] for box in boxes if box[1] == col), default=1)
        for col in range(n_cols)
    ]
    slot_h = [
        max((box[2] + box[4] for box in boxes if box[0] == row), default=1)
        for row in range(n_rows)
    ]
    bar_frac = 0.11  # colorbar width, as a fraction of the slot left of it
    gap_frac = 0.09  # white space between neighbouring slots in one row
    col_in = [w * cell_width / dpi for w in slot_w]
    bar_in = [bar_frac * w for w in col_in]
    gap_in = [gap_frac * w for w in col_in]
    band_in = [h * cell_height / dpi for h in slot_h]
    width_in = left_in + sum(col_in) + sum(bar_in) + sum(gap_in[: n_cols - 1])
    height_in = top_in + bottom_in + sum(band_in)
    # the constrained layout gives every label the margin it needs, which the
    # hand-placed rectangles alone would often cut off; to take part in that
    # layout, the panels and their colorbars are subplots of a grid whose
    # column widths and row heights keep the pixel proportions of the figure
    figure = plt.figure(figsize=(width_in, height_in), layout="constrained")
    grid = gridspec.GridSpec(
        n_rows,
        2 * n_cols,
        figure=figure,
        width_ratios=[
            width for col in range(n_cols) for width in (col_in[col], bar_in[col])
        ],
        height_ratios=band_in,
        wspace=gap_frac,
        hspace=0.02,  # the thin white margin between neighbouring row bands
    )

    # the panel axes keep the last row after the drawing loop below; those are
    # the ones naming the columns of the whole figure
    panel_axes: dict[int, object] = {}
    # the leftmost panel axes of every row band, to name that band
    left_axes: dict[int, object] = {}
    # draw the panels row band by row band; the cells of a band are colored by
    # its colorbar, stretched onto the data values of its two ends, while every
    # panel still gets a colorbar of its own, as the recreation does not know
    # that the bars of a band are shared
    for row in range(n_rows):
        row_values = np.asarray(bars[row]["values"], dtype=float)
        vmin, vmax = float(row_values.min()), float(row_values.max())
        # the extracted colorbar runs from its top end to its bottom end, while
        # a colormap runs from its smallest value upwards, so a colorbar whose
        # top carries the largest value has to be reversed to color the panels
        descending = row_values[0] > row_values[1]
        colors = colorlines[row][::-1] if descending else colorlines[row]
        cmap = ListedColormap(colors / 255)
        cmap.set_bad("white")
        mappable = plt.cm.ScalarMappable(
            cmap=cmap, norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        )
        for col in range(n_cols):
            axes = figure.add_subplot(grid[row, 2 * col])
            panel_axes[col] = axes
            if col == 0:
                left_axes[row] = axes
            # the panels carry no ticks or frame, but their axes stay visible so
            # that the labels naming them take part in the constrained layout
            axes.set_xticks([])
            axes.set_yticks([])
            axes.patch.set_visible(False)
            axes.spines[:].set_visible(False)
            for r, c, oy, ox, ny, nx in boxes:
                if (r, c) != (row, col) or ny <= 0 or nx <= 0:
                    continue
                axes.imshow(
                    stack[oy : oy + ny, ox : ox + nx, row, col],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    extent=[ox, ox + nx, oy + ny, oy],
                    aspect="auto",
                    interpolation="nearest",
                )
            axes.set_xlim(0, slot_w[col])
            axes.set_ylim(slot_h[row], 0)
            bar_axes = figure.add_subplot(grid[row, 2 * col + 1])
            colorbar = figure.colorbar(mappable, cax=bar_axes)
            colorbar.ax.tick_params(labelsize=6, pad=1.5, length=1.5)
    # the names of the rows, columns, and the figure, drawn as the labels of
    # the axes they name so that the constrained layout keeps them readable
    for row in range(n_rows):
        if row_labels and row_labels[row] and row in left_axes:
            left_axes[row].set_ylabel(
                row_labels[row],
                rotation=0,
                va="center",
                ha="right",
                fontsize="small",
                labelpad=4,
            )
    if column_labels:
        for col in range(n_cols):
            if col in panel_axes:
                panel_axes[col].set_xlabel(column_labels[col], fontsize="small")
    if title:
        figure.suptitle(title)
    figure.savefig(filename, dpi=dpi * save_scale)
    return figure


def main() -> None:
    """Digitize every panel of the figure and stack the values."""
    prs = argparse.ArgumentParser(description=__doc__ or "")
    prs.add_argument(
        "lattice_json", type=str, help="Lattice file written by panel_lattice.py."
    )
    prs.add_argument(
        "bars_json",
        type=str,
        help=(
            "JSON file with one 'bars' entry per panel row: a legend_box "
            "[left, top, right, bottom] with left = right at the colorbar of "
            "that row, and the data values [top, bottom] at its two ends."
        ),
    )
    prs.add_argument(
        "output",
        type=str,
        help="Filename of the stacked data array, with a .npy extension.",
    )
    prs.add_argument(
        "--save-panels",
        type=str,
        metavar="DIRECTORY",
        default=None,
        help="Also save every panel's values as DIRECTORY/<stem>_r<row>_c<col>.txt.",
    )
    prs.add_argument(
        "--recreate-image",
        type=str,
        metavar="FILENAME",
        default=None,
        help=(
            "If provided, also save the recreated figure to this filename, with "
            "a small colorbar next to every panel as in the original figure."
        ),
    )
    prs.add_argument(
        "--row-labels",
        type=str,
        default=None,
        metavar="JSON",
        help=(
            "JSON-formatted string of the label of every panel row, e.g. "
            '["SO2 [ppmv]", ...]; drawn to the left of each row of panels.'
        ),
    )
    prs.add_argument(
        "--column-labels",
        type=str,
        default=None,
        metavar="JSON",
        help=(
            "JSON-formatted string of the label of every panel column, e.g. "
            '["a", "b", "c", "d"]; drawn below each column of panels.'
        ),
    )
    args = prs.parse_args()
    assert args.output.endswith(".npy"), "the stacked array needs a .npy filename"

    with open(args.lattice_json, encoding="utf-8") as handle:
        lattice = json.load(handle)
    with open(args.bars_json, encoding="utf-8") as handle:
        bars = json.load(handle)["bars"]
    n_rows, n_cols = lattice["n_rows"], lattice["n_cols"]
    assert len(bars) == n_rows, f"expected {n_rows} colorbar entries, got {len(bars)}"
    for bar in bars:
        check_bar(bar)
    row_labels = (
        [str(label) for label in json.loads(args.row_labels)]
        if args.row_labels is not None
        else None
    )
    if row_labels is None:
        # the rows of a figure are usually the variables of its colorbars,
        # which the bars file already names
        named = [str(bar.get("label", "")) for bar in bars]
        row_labels = named if any(named) else None
    column_labels = (
        [str(label) for label in json.loads(args.column_labels)]
        if args.column_labels is not None
        else None
    )
    assert (
        row_labels is None or len(row_labels) == n_rows
    ), f"expected {n_rows} row labels, got {len(row_labels)}"
    assert (
        column_labels is None or len(column_labels) == n_cols
    ), f"expected {n_cols} column labels, got {len(column_labels)}"

    img = load_image(lattice["image"]).astype(np.int64)
    cell_width = lattice["cell_size_px"]["width"]
    cell_height = lattice["cell_size_px"]["height"]
    panel_h = int(round(lattice["panel_height_px"]))
    grids = {
        (panel["row"], panel["col"]): panel for panel in lattice["cell_grid"]["panels"]
    }

    # absolute seam borders of every panel, extended over its window
    borders: dict[tuple[int, int], tuple[NDArray[np.int64], NDArray[np.int64]]] = {}
    for panel in lattice["panels"]:
        key = (panel["row"], panel["col"])
        panel_grid = grids[key]
        window_x0 = int(panel_grid["x0"])
        window_y0 = int(panel_grid["y0"])
        borders[key] = (
            window_borders(
                panel_grid,
                1,
                window_x0,
                window_x0 + int(panel_grid["width"]),
                cell_width,
            ),
            window_borders(
                panel_grid,
                0,
                window_y0,
                min(window_y0 + panel_h, img.shape[0]),
                cell_height,
            ),
        )

    # the earliest seam of every column and row band is the origin of the
    # stacked indices of that column or band; the stack extends as far as the
    # data of the widest column or row band reaches
    origin_x = {
        col: min(borders[(row, col)][0][0] for row in range(n_rows))
        for col in range(n_cols)
    }
    origin_y = {
        row: min(borders[(row, col)][1][0] for col in range(n_cols))
        for row in range(n_rows)
    }
    nx_big = max(
        int(np.ceil((borders[key][0][-1] - origin_x[key[1]]) / cell_width))
        for key in borders
    )
    ny_big = max(
        int(np.ceil((borders[key][1][-1] - origin_y[key[0]]) / cell_height))
        for key in borders
    )
    values = np.full((ny_big, nx_big, n_rows, n_cols), np.nan)
    print(f"stack shape (y, x, row, column): {values.shape}")

    masked_total = 0
    boxes: list[tuple[int, int, int, int, int, int]] = []
    for row in range(n_rows):
        for col in range(n_cols):
            key = (row, col)
            borders_x, borders_y = borders[key]
            offset_x, offset_y, nx, ny = panel_placement(
                borders_x,
                borders_y,
                origin_x[col],
                origin_y[row],
                cell_width,
                cell_height,
                nx_big,
                ny_big,
            )
            panel_values, n_masked = extract_panel(img, bars[row], borders_x, borders_y)
            values[offset_y : offset_y + ny, offset_x : offset_x + nx, row, col] = (
                panel_values[:ny, :nx]
            )
            boxes.append((row, col, offset_y, offset_x, ny, nx))
            masked_total += n_masked
            print(
                f"panel ({row}, {col}): {ny} x {nx} cells into slot "
                f"(y {offset_y}, x {offset_x}), {n_masked} masked cells"
            )
            if args.save_panels is not None:
                os.makedirs(args.save_panels, exist_ok=True)
                stem = args.lattice_json.rpartition(".")[0]
                save_values(
                    panel_values,
                    os.path.join(args.save_panels, f"{stem}_r{row}_c{col}.txt"),
                    False,
                )

    np.save(args.output, values)
    n_data = int(np.count_nonzero(~np.isnan(values)))
    print(
        f"wrote {args.output}: {n_data} values, "
        f"{values.size - n_data - masked_total} empty slots, {masked_total} "
        f"masked cells"
    )

    if args.recreate_image is not None:
        recreate_figure(
            stack=values,
            boxes=boxes,
            bars=bars,
            colorlines=[colorline_of_bar(img, bar["legend_box"]) for bar in bars],
            cell_width=cell_width,
            cell_height=cell_height,
            row_labels=row_labels,
            column_labels=column_labels,
            title=os.path.basename(lattice["image"]),
            filename=args.recreate_image,
            dpi=150.0,
            save_scale=2.0,
        )
        print(f"wrote {args.recreate_image}")


if __name__ == "__main__":
    main()
