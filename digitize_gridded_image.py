"""
Given an image file, an outcrop, and a legend area, convert the colors
in a regular gridded mesh back to the data space.
"""

import sys
import json
from collections.abc import Callable
from typing import Any, Literal

import numpy as np
from argparse import ArgumentParser
from numpy.typing import NDArray
from PIL import Image
from scipy.spatial.distance import cdist, euclidean
from scipy.interpolate import CubicSpline

# An RGB image as an array of shape (height, width, 3) with integer color
# channels in the range 0-255.
RGBImage = NDArray[Any]
# Boxes of inclusive pixel boundary coordinates in the order
# [left, top, right, bottom].
Box = NDArray[Any] | tuple[int, int, int, int]
# Extracted data values of a grid of cells; either floats for a continuous
# colorbar or categorical labels such as strings.
Values = NDArray[Any]
# Maps a position along a colorbar (as a pixel-like index) to a data value.
IndexToValue = Callable[[NDArray[Any]], Values]


# copied from https://stackoverflow.com/a/30305181/13054094,
# licensed under https://creativecommons.org/licenses/by-sa/3.0/,
# below the original comment from orlp:
# ------------------------------------------------------------------------------
# I implemented Yehuda Vardi and Cun-Hui Zhang's algorithm for the geometric
# median, described in their paper "The multivariate L1-median and associated
# data depth". Everything is vectorized in numpy, so should be very fast. I
# didn't implement weights - only unweighted points.
# In addition to the default SO license terms, I release the code above under
# the zlib license, if you so prefer.
# ------------------------------------------------------------------------------
def geometric_median(X, eps=1e-5):
    y = np.mean(X, 0)

    while True:
        D = cdist(X, [y])
        nonzeros = (D != 0)[:, 0]

        Dinv = 1 / D[nonzeros]
        Dinvs = np.sum(Dinv)
        W = Dinv / Dinvs
        T = np.sum(W * X[nonzeros], 0)

        num_zeros = len(X) - np.sum(nonzeros)
        if num_zeros == 0:
            y1 = T
        elif num_zeros == len(X):
            return y
        else:
            R = (T - y) * Dinvs
            r = np.linalg.norm(R)
            rinv = 0 if r == 0 else num_zeros / r
            y1 = max(0, 1 - rinv) * T + min(1, rinv) * y

        if euclidean(y, y1) < eps:
            return y1

        y = y1


def project_point_on_line_3d(
    line: np.ndarray, point: np.ndarray
) -> tuple[int, np.ndarray]:
    """
    Find the vertex on a line made from multiple segments that is closest to another
    point.
    Based on https://math.stackexchange.com/q/1521128 with further modifications to
    accept a segmented line and to restrict the output to be one of the input vertices.

    Parameters
    ----------
    line
        3D vertices that make up the segmented line.
    point
        Point to find the closest vertex to.

    Returns
    -------
    index
        The index of the input vertex closest to the input point.
    vertex
        The coordinates of the closest vertex.
    """
    # input checks
    # line needs to be 2D, each row are 3D vertex coordinates
    # the point needs to be 3D
    assert line.ndim == 2
    assert point.ndim == 1
    assert line.shape[0] > 1
    assert line.shape[1] == point.size == 3
    # reshape for readability later
    p = point.reshape(1, 3)
    # get the Q (starting) and R (ending) points of individual segments
    q = np.atleast_2d(line[:-1, :])
    r = np.atleast_2d(line[1:, :])
    # get segment vectors
    r_minus_q = r - q
    # get vectors to point
    q_minus_p = q - p
    # initialize loop over line segments
    t = np.empty(line.shape[0] - 1)
    g = np.empty((t.size, 3))
    for i in range(t.size):
        # get distance from starting point to point on line that is closest
        # to the target point if the line segment were infiniely long
        that = np.dot(r_minus_q[i, :], q_minus_p[i, :]) / np.dot(
            r_minus_q[i, :], r_minus_q[i, :]
        )
        # restrict that distance to yield only points on the segment itself
        t[i] = min(max(that, 0), -1)
        # save the point on the segment that is closest to the target point
        g[i, :] = q[i, :] - t[i] * r_minus_q[i, :]
    # the globally closest point on the line has the minimum Euclidean distance
    dists_line = np.sqrt(np.sum((g - p) ** 2, axis=1))
    # find index and point closest to target point on entire line
    ix_line = np.argmin(dists_line)
    g_line = g[ix_line, :]
    # if it's not already the vertex, find the vertex that is closest to the point
    # on the line segment
    dists_vertices = np.sqrt(np.sum((g_line - line) ** 2, axis=1))
    # return the index and location of the vertex that is closest to the target point
    index = np.argmin(dists_vertices)
    vertex = line[index, :]
    return index, vertex


def load_image(filename: str) -> RGBImage:
    """
    Read an image file as an array of integer RGB colors.

    Parameters
    ----------
    filename
        Filename of the image to read.

    Returns
    -------
    img
        Array of shape (height, width, 3) holding the integer RGB channels of
        the image, with any opacity values discarded.
    """
    return np.asarray(Image.open(filename))[:, :, :3].astype(int)


def cell_borders(box: Box, n: int, axis: Literal["x", "y"]) -> NDArray[Any]:
    """
    Evenly divide a bounding box into a regular grid of cells along one axis.

    Parameters
    ----------
    box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        area to divide.
    n
        Number of cells along that axis.
    axis
        Either "x" to divide along the horizontal or "y" along the vertical
        direction of the box.

    Returns
    -------
    borders
        The `n` + 1 integer pixel coordinates delimiting the cells, where the
        last border lies one pixel past the edge of the bounding box so that
        it can be used as the exclusive end of a slice.
    """
    lo, hi = (box[0], box[2]) if axis == "x" else (box[1], box[3])
    return np.linspace(lo, hi + 1, num=n + 1).round().astype(int)


def extract_dominant_colors(
    img: RGBImage,
    data_box: Box,
    nx: int,
    ny: int,
    ignore_black: int | None = None,
) -> NDArray[Any]:
    """
    Extract the dominant color of every cell of a regular grid of cells.

    Each cell color is estimated with something similar to the median but in
    three dimensions, which is robust against the mixed colors occurring where
    neighboring cells meet.

    Parameters
    ----------
    img
        Image array with integer RGB channels.
    data_box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        data field, i.e. of the grid of colored cells.
    nx
        Number of data grid points in the x direction.
    ny
        Number of data grid points in the y direction.
    ignore_black
        If provided, pixels this far away in RGB space from true black (using
        Euclidean 3D distances) are excluded from the dominant color of a cell,
        which serves to ignore masked-out or background pixels.

    Returns
    -------
    cmed
        Array of shape (ny, nx, 3) with the dominant integer RGB color of
        every cell.
    """
    # get grid that defines the different sub-images
    borders_x = cell_borders(data_box, nx, "x")
    borders_y = cell_borders(data_box, ny, "y")
    # initialize loop over sub-images
    cmed = np.empty((borders_y.size - 1, borders_x.size - 1, 3))
    for i in range(cmed.shape[0]):
        for j in range(cmed.shape[1]):
            # extract sub-image
            subimg = img[
                borders_y[i] : borders_y[i + 1], borders_x[j] : borders_x[j + 1], :
            ].reshape(-1, 3)
            # if ignoring black colors, remove the pixels that are too close
            # given the input radius
            if ignore_black is not None:
                not_black_ix = (
                    np.sqrt(np.sum(subimg.astype(float) ** 2, axis=1)) > ignore_black
                )
                subimg = subimg[not_black_ix, :]
            # get the most common ("dominant") color using something similar to
            # the median but in 3D
            cmed[i, j, :] = geometric_median(subimg)
    # make sure the colors are in the normal 0-255 integer RGB range
    return np.clip(cmed, a_min=0, a_max=255).astype(int)


def extract_colorline(img: RGBImage, legend_box: Box) -> NDArray[Any]:
    """
    Extract the colorbar of an image as a sequence of colors.

    Parameters
    ----------
    img
        Image array with integer RGB channels.
    legend_box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        colorbar, which has to be either perfectly vertical or horizontal.

    Returns
    -------
    colorline
        Array of shape (length, 3) with the integer RGB color of every pixel
        along the colorbar, ordered from top to bottom if the colorbar is
        vertical and from left to right if it is horizontal.
    """
    return img[
        legend_box[1] : legend_box[3] + 1, legend_box[0] : legend_box[2] + 1, :
    ].squeeze()


def marker_colorline_indices(
    markers_pixel: NDArray[Any], legend_box: Box, n_colorline: int
) -> NDArray[Any]:
    """
    Map the pixel coordinates of the colorbar markers onto the colorbar.

    Parameters
    ----------
    markers_pixel
        Ascending pixel coordinates of the tick labels along the colorbar.
    legend_box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        colorbar.
    n_colorline
        Number of pixels that make up the extracted colorbar.

    Returns
    -------
    marker_indices
        The index into the extracted colorbar of every marker.
    """
    # map the pixel coordinates linearly onto the indices of the colorbar
    return (
        (
            (markers_pixel - legend_box[1])
            / (legend_box[3] - legend_box[1])
            * (n_colorline - 1)
        )
        .round()
        .astype(int)
    )


def unique_colors(colorline: NDArray[Any]) -> NDArray[Any]:
    """
    Find the distinct colors of a colorbar, preserving their order.

    Parameters
    ----------
    colorline
        Array of shape (length, 3) with the colors along a colorbar.

    Returns
    -------
    colorline_unique
        Array of shape (n_unique, 3) with the colors that appear on the
        colorbar, in the order in which they first appear.
    """
    # the first occurrence of every color defines its position in the sequence
    _, colorline_uniqind = np.unique(colorline, axis=0, return_index=True)
    return colorline[np.sort(colorline_uniqind)]


def colors_to_colorline_indices(
    cmed: NDArray[Any], colorline: NDArray[Any], colorline_unique: NDArray[Any]
) -> tuple[NDArray[Any], NDArray[Any]]:
    """
    Locate the dominant color of every cell of a grid on the colorbar.

    Parameters
    ----------
    cmed
        Array of shape (ny, nx, 3) with the dominant color of every cell, as
        returned by `extract_dominant_colors`.
    colorline
        Array of shape (length, 3) with the colors along the colorbar.
    colorline_unique
        Array of shape (n_unique, 3) with the distinct colors along the
        colorbar, as returned by `unique_colors`.

    Returns
    -------
    indices_unique
        Array of shape (ny, nx) with the index into the distinct colors of the
        colorbar that lies closest to the dominant color of every cell.
    indices_general
        Array of shape (ny, nx) with the mean index into the pixels of the
        colorbar whose color matches the dominant color of every cell.
    """
    # initialize loop over all dominant colors
    indices_unique = np.empty(cmed.shape[:2], dtype=int)
    indices_general = np.empty_like(indices_unique)
    for i in range(cmed.shape[0]):
        for j in range(cmed.shape[1]):
            # find the color on the colorbar that is closest in 3D RGB space
            # to the dominant color that was extracted in the previous step
            indices_unique[i, j], gix = project_point_on_line_3d(
                colorline_unique, cmed[i, j, :]
            )
            # since colors can repeat in a colorbar, we need to decide which pixel
            # to map into data space - here we use the mean one
            indices_general[i, j] = np.flatnonzero(
                np.all(gix == colorline, axis=1)
            ).mean()
    return indices_unique, indices_general


def sorted_direction(markers_data: NDArray[Any]) -> int:
    """
    Determine whether the data values of the markers ascend or descend.

    Parameters
    ----------
    markers_data
        Data values of the tick labels of a colorbar.

    Returns
    -------
    direction
        1 if the values are ascending and -1 if they are descending.

    Raises
    ------
    AssertionError
        If the values are neither ascending nor descending.
    """
    if np.all(np.sort(markers_data) == markers_data):
        return 1
    if np.all(np.sort(markers_data) == np.flip(markers_data)):
        return -1
    raise AssertionError


def build_index_to_value(
    marker_indices: NDArray[Any], markers_data: NDArray[Any]
) -> tuple[IndexToValue, bool]:
    """
    Build the function that maps positions along a colorbar to data values.

    Parameters
    ----------
    marker_indices
        Index into the colorbar of every marker, as returned by
        `marker_colorline_indices`.
    markers_data
        Data values of the markers, corresponding to `marker_indices`.

    Returns
    -------
    index_to_value
        Function that maps an array of colorbar indices to data values. It
        expects a two-dimensional array of the shape of the grid of cells and
        returns values of that same shape.
    categorical
        False if the markers describe a continuous colorbar, in which case the
        mapping interpolates cubically between the markers. True if the markers
        are categorical labels, in which case the mapping returns the label of
        the closest marker for a flat array of colorbar indices.
    """
    try:
        return CubicSpline(marker_indices, markers_data), False
    # catch the case where the markers are categorical
    except ValueError:
        pass

    def index_to_value(index: NDArray[Any]) -> Values:
        # give every colorbar index the value of its closest marker
        closest = np.argmin(
            np.abs(marker_indices[:, None] - index.reshape(1, -1)), axis=0
        )
        return markers_data[closest]

    return index_to_value, True


def map_indices_to_values(
    index_to_value: IndexToValue,
    categorical: bool,
    indices_general: NDArray[Any],
    shape: tuple[int, int],
) -> Values:
    """
    Convert the colorbar locations of a grid of cells to data values.

    Parameters
    ----------
    index_to_value
        Function that maps colorbar indices to data values, as returned by
        `build_index_to_value`.
    categorical
        Whether that mapping describes categorical labels rather than a
        continuous colorbar.
    indices_general
        Array of shape (ny, nx) with the colorbar location of every cell, as
        returned by `colors_to_colorline_indices`.
    shape
        The shape (ny, nx) of the grid of cells, which the lookup of categorical
        labels needs in order to reshape its output.

    Returns
    -------
    values
        Array of shape (ny, nx) with the data value of every cell.
    """
    if categorical:
        # the lookup of labels works on a flat array and has to be reshaped
        return index_to_value(indices_general).reshape(shape)
    return index_to_value(indices_general)


def save_values(values: Values, output_name: str, categorical: bool) -> None:
    """
    Save a grid of extracted data values.

    Parameters
    ----------
    values
        The data values to save.
    output_name
        Output filename with either a `.npy` or a `.txt` extension.
    categorical
        Whether the values are categorical labels rather than floats, which
        decides the text format used for `.txt` output.

    Raises
    ------
    AssertionError
        If the filename has neither a `.npy` nor a `.txt` extension.
    """
    assert ("." in output_name) and (output_name.split(".")[-1] in ["npy", "txt"])
    out_fmt = output_name.split(".")[-1]
    if out_fmt == "npy":
        np.save(output_name, values)
    elif out_fmt == "txt":
        np.savetxt(output_name, values, fmt="%s" if categorical else "%.18e")


def recreate_image(
    values: Values,
    categorical: bool,
    indices_general: NDArray[Any],
    colorline: NDArray[Any],
    marker_indices: NDArray[Any],
    markers_data: NDArray[Any],
    filename: str,
    marker_data_sorted: int | None = None,
) -> None:
    """
    Save a debugging image that recreates a grid of cells from its colors.

    Parameters
    ----------
    values
        The extracted data values of every cell.
    categorical
        Whether those values are categorical labels rather than floats.
    indices_general
        Array of shape (ny, nx) with the colorbar location of every cell, as
        returned by `colors_to_colorline_indices`.
    colorline
        Array of shape (length, 3) with the colors along the colorbar.
    marker_indices
        Index into the colorbar of every marker, as returned by
        `marker_colorline_indices`.
    markers_data
        Data values of the markers, corresponding to `marker_indices`.
    filename
        Filename under which to save the recreated image.
    marker_data_sorted
        1 if `markers_data` ascends and -1 if it descends. Defaults to None, in
        which case the direction is determined from `markers_data`.
    """
    # delayed imports
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.cm import ScalarMappable

    if marker_data_sorted is None:
        marker_data_sorted = sorted_direction(markers_data)

    # create function that maps data to original colorbar indices
    if categorical:
        # reuse the step indices of our markers directly
        # make a discrete colormap using those values
        unique_colindex = marker_indices
        borders = np.array(
            [0]
            + ((unique_colindex[1:] + unique_colindex[:-1]) / 2).tolist()
            + [colorline.shape[0] - 1]
        )
    else:
        if marker_data_sorted == 1:  # ascending
            value_to_index = CubicSpline(markers_data, marker_indices)
        else:  # descending
            value_to_index = CubicSpline(np.flip(markers_data), np.flip(marker_indices))
        # get the unique values and colors we extracted
        unique_values = np.unique(values.ravel())
        unique_colindex = value_to_index(unique_values).round().astype(int)
        # make a discrete colormap using those values
        if marker_data_sorted == 1:
            borders = np.array(
                [markers_data[0]]
                + ((unique_values[1:] + unique_values[:-1]) / 2).tolist()
                + [markers_data[-1]]
            )
        else:
            borders = np.array(
                [markers_data[-1]]
                + ((unique_values[1:] + unique_values[:-1]) / 2).tolist()
                + [markers_data[0]]
            )
    colors_unique = colorline[unique_colindex]
    cmap = ListedColormap(colors_unique / 255)
    norm = BoundaryNorm(borders, cmap.N)
    # start the figure
    plt.figure()
    # add the recreated data using the colors we estimated
    plt.imshow(
        indices_general.astype(float) if categorical else values,
        cmap=cmap,
        norm=norm,
    )
    # add the colorbar
    cbar = plt.colorbar(
        mappable=ScalarMappable(norm=norm, cmap=cmap),
        ax=plt.gca(),
        label="Data Units",
    )
    if categorical:
        cbar.set_ticks(marker_indices.astype(float), labels=markers_data)
        cbar.ax.invert_yaxis()
    else:
        cbar.set_ticks(unique_values, labels=[""] * unique_values.size, minor=True)
        cbar.set_ticks(
            unique_values[np.linspace(0, unique_values.size - 1, num=10, dtype=int)]
        )
    # add info
    plt.title("Recreated Image")
    plt.xlabel("X Index")
    plt.ylabel("Y Index")
    # save
    plt.savefig(filename)


def validate_inputs(
    data_box: Box,
    legend_box: Box,
    markers_pixel: NDArray[Any],
    markers_data: NDArray[Any],
    output_name: str,
    ignore_black: int | None,
) -> int:
    """
    Check the inputs of the extraction and report the direction of the colorbar.

    Parameters
    ----------
    data_box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        data field.
    legend_box
        Inclusive pixel boundary coordinates [left, top, right, bottom] of the
        colorbar.
    markers_pixel
        Ascending pixel coordinates of the tick labels along the colorbar.
    markers_data
        Data values of the tick labels.
    output_name
        Output filename for data with either a `.npy` or a `.txt` extension.
    ignore_black
        Distance in RGB space up to which pixels are treated as black, or None
        if that filter is unset.

    Returns
    -------
    marker_data_sorted
        1 if the data values of the markers ascend and -1 if they descend.

    Raises
    ------
    AssertionError
        If any of the inputs is inconsistent.
    """
    # each bounding box has to have four elements that are ordered ascending
    assert data_box.size == 4
    assert (data_box[2] > data_box[0]) and (data_box[3] > data_box[1])
    assert legend_box.size == 4
    assert (legend_box[2] >= legend_box[0]) and (legend_box[3] >= legend_box[1])
    # the legend bounding box has to be either perfectly vertical or horizontal
    assert np.logical_xor(
        legend_box[0] == legend_box[2], legend_box[1] == legend_box[3]
    )
    # there have to be the same number of pixel and value markers for the colorbar
    assert markers_pixel.size == markers_data.size
    # the markers need to be sorted
    assert np.all(np.sort(markers_pixel) == markers_pixel)
    # the data needs to be sorted but the direction doesn't matter
    marker_data_sorted = sorted_direction(markers_data)
    # test for a valid filename with extension
    assert ("." in output_name) and (output_name.split(".")[-1] in ["npy", "txt"])
    # test that ignore_black is either unset or non-negative
    assert (ignore_black is None) or (ignore_black >= 0)
    return marker_data_sorted


def main() -> None:
    """
    Run the extraction of the data values of a gridded image from its colorbar.
    """
    # make argument parser
    prs = ArgumentParser(description="Extract data from an image and colorbar.")
    prs.add_argument("filename", type=str, help="Filename of image.")
    prs.add_argument(
        "data_box",
        type=str,
        help=(
            "JSON-formatted string of (inclusive) pixel boundary coordinates "
            "of the data in the order [left, top, right, bottom]."
        ),
    )
    prs.add_argument(
        "legend_box",
        type=str,
        help=(
            "JSON-formatted string of (inclusive) pixel boundary coordinates "
            "of the colorbar in the order [left, top, right, bottom]; "
            "note that either left=right or top=bottom."
        ),
    )
    prs.add_argument(
        "markers_pixel",
        type=str,
        help=(
            "JSON-formatted string of ascending pixel coordinates of the "
            "tick labels (x coordinates if colorbar horzontal, y otherwise)."
        ),
    )
    prs.add_argument(
        "markers_data",
        type=str,
        help="JSON-formatted string of data values of the tick labels.",
    )
    prs.add_argument(
        "nx", type=int, help="Number of data grid points in the x direction."
    )
    prs.add_argument(
        "ny",
        type=int,
        help="Number of data grid points in the y direction (defaults to `nx`).",
    )
    prs.add_argument(
        "output_name",
        type=str,
        help="Output filename for data with either `.npy` or `.txt` extension.",
    )
    prs.add_argument(
        "--recreate-image",
        type=str,
        default=None,
        help=(
            "If provided, also save a debugging image with the extracted colors "
            "to the specified filename."
        ),
    )
    prs.add_argument(
        "--ignore-black",
        type=int,
        default=None,
        help=(
            "When extracting the dominant colors, ignore pixels this far away in "
            "RGB space from true black using Euclidean 3D distances."
        ),
    )

    # print help if wrongly called
    if len(sys.argv) == 1:
        prs.print_help()
        exit()

    # parse arguments
    args = prs.parse_args()
    nx = args.nx
    ny = args.ny
    data_box = np.asarray(json.loads(args.data_box)).ravel()
    legend_box = np.asarray(json.loads(args.legend_box)).ravel()
    markers_pixel = np.asarray(json.loads(args.markers_pixel)).ravel()
    markers_data = np.asarray(json.loads(args.markers_data)).ravel()
    output_name = args.output_name
    recreate = args.recreate_image
    ignore_black = args.ignore_black

    # input checks
    marker_data_sorted = validate_inputs(
        data_box, legend_box, markers_pixel, markers_data, output_name, ignore_black
    )

    # load image file as NumPy array, discard opacity values
    img = load_image(args.filename)

    # get median colors from each grid value
    cmed = extract_dominant_colors(img, data_box, nx, ny, ignore_black)

    # get best-fit data value from colors
    colorline = extract_colorline(img, legend_box)
    # map the pixel indices of the markers to the colorbar location
    marker_indices = marker_colorline_indices(
        markers_pixel, legend_box, colorline.shape[0]
    )
    # find the unique colors of the colorbar and preserve their order
    colorline_unique = unique_colors(colorline)
    indices_unique, indices_general = colors_to_colorline_indices(
        cmed, colorline, colorline_unique
    )
    # create an interpolant function that maps colorbar indices to the data space
    index_to_value, categorical = build_index_to_value(marker_indices, markers_data)
    # convert the extracted colorbar indices to data space
    values = map_indices_to_values(
        index_to_value, categorical, indices_general, (ny, nx)
    )

    # save 2D data array
    save_values(values, output_name, categorical)

    # recreate image if desired
    if recreate is not None:
        recreate_image(
            values=values,
            categorical=categorical,
            indices_general=indices_general,
            colorline=colorline,
            marker_indices=marker_indices,
            markers_data=markers_data,
            filename=recreate,
            marker_data_sorted=marker_data_sorted,
        )


if __name__ == "__main__":
    main()
