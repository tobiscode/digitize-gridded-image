"""
Given an image file, an outcrop, and a legend area, convert the colors
in a regular gridded mesh back to the data space.
"""

import sys
import json
import numpy as np
from argparse import ArgumentParser
from PIL import Image
from scipy.spatial.distance import cdist, euclidean
from scipy.interpolate import CubicSpline


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


if __name__ == "__main__":

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
    recreate_image = args.recreate_image
    ignore_black = args.ignore_black

    # input checks
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
    if np.all(np.sort(markers_data) == markers_data):
        marker_data_sorted = 1
    elif np.all(np.sort(markers_data) == np.flip(markers_data)):
        marker_data_sorted = -1
    else:
        raise AssertionError
    # test for a valid filename with extension
    assert ("." in output_name) and (output_name.split(".")[-1] in ["npy", "txt"])
    # test that ignore_black is either unset or non-negative
    assert (ignore_black is None) or (ignore_black >= 0)

    # load image file as NumPy array, discard opacity values
    img = np.asarray(Image.open(args.filename))[:, :, :3].astype(int)

    # get median colors from each grid value
    # get grid that defines the different sub-images
    borders_x = (
        np.linspace(data_box[0], data_box[2] + 1, num=nx + 1).round().astype(int)
    )
    borders_y = (
        np.linspace(data_box[1], data_box[3] + 1, num=ny + 1).round().astype(int)
    )
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
    cmed = np.clip(cmed, a_min=0, a_max=255).astype(int)

    # get best-fit data value from colors
    # get colorbar
    colorline = img[
        legend_box[1] : legend_box[3] + 1, legend_box[0] : legend_box[2] + 1, :
    ].squeeze()
    # map the pixel indices of the markers to the colorbar location
    steps_plot_indicized = (
        (
            (markers_pixel - legend_box[1])
            / (legend_box[3] - legend_box[1])
            * (colorline.shape[0] - 1)
        )
        .round()
        .astype(int)
    )
    # find the unique colors of the colorbar and preserve their order
    _, colorline_uniqind = np.unique(colorline, axis=0, return_index=True)
    colorline_unique = colorline[np.sort(colorline_uniqind)]
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
    # create an interpolant function that maps colorbar indices to the data space
    try:
        pixel_to_value = CubicSpline(steps_plot_indicized, markers_data)
        categorical = False
    # catch the case where the markers are categorical
    except ValueError:
        pixel_to_value = lambda index: markers_data[
            np.argmin(
                np.abs(steps_plot_indicized[:, None] - index.reshape(1, -1)), axis=0
            ).reshape(args.ny, args.nx)
        ]
        categorical = True
    # convert the extracted colorbar indices to data space
    values = pixel_to_value(indices_general)

    # save 2D data array
    out_fmt = output_name.split(".")[-1]
    if out_fmt == "npy":
        np.save(output_name, values)
    elif out_fmt == "txt":
        np.savetxt(output_name, values, fmt="%s" if categorical else "%.18e")

    # recreate image if desired
    if recreate_image is not None:
        # delayed imports
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap, BoundaryNorm
        from matplotlib.cm import ScalarMappable

        # create function that maps data to original colorbar indices
        if categorical:
            # reuse the step indices of our markers directly
            # make a discrete colormap using those values
            unique_colindex = steps_plot_indicized
            borders = np.array(
                [0]
                + ((unique_colindex[1:] + unique_colindex[:-1]) / 2).tolist()
                + [colorline.shape[0] - 1]
            )
        else:
            if marker_data_sorted == 1:  # ascending
                value_to_index = CubicSpline(markers_data, steps_plot_indicized)
            else:  # descending
                value_to_index = CubicSpline(
                    np.flip(markers_data), np.flip(steps_plot_indicized)
                )
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
        unique_colors = colorline[unique_colindex]
        cmap = ListedColormap(unique_colors / 255)
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
            cbar.set_ticks(steps_plot_indicized.astype(float), labels=markers_data)
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
        plt.savefig(recreate_image)
