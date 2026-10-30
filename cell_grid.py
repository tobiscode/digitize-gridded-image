"""
Fit the integer cell grid inside every panel of a multi-panel figure.

The width and height of the cells are shared by all panels and already
measured by ``cell_fft.py`` (stored in the lattice JSON). Only the phase of
the seam lattice inside each panel is left open. Starting from the panel
window bounding boxes of the lattice, the nominal seam positions inside a
window are its multiples of the shared cell size, shifted by a sub-pixel
offset: the column windows
follow the white gutters of the panel content and so hold no fixed anchor
offset, while the row windows do and pass that offset in as the nominal start.
A seam is an antialiased boundary between two cells, so the pixels on both
sides of it differ in colour while cell interiors are flat; a candidate seam
lattice is therefore scored by the MEAN absolute RGB difference across its
seam lines (mean rather than sum, so that a phase which rounds one extra
fringe line into the window gains nothing), and white or masked areas simply
contribute nothing. The script probes sub-pixel shifts of the nominal seams
(less than one cell, followed by rounding to integer pixel seams), picks the
best per panel and per axis, and reports for each choice the strongest
competing phase of the score curve together with the fraction of the chosen
seam lines that are individually strong, which exposes a phase that only one
strong stray structure supports. It also reports how strongly the panel phases
of a figure agree as one global lattice, stores the integer seam positions in
the JSON, and can draw the fitted cell grid over the figure as a debug image.
"""

import argparse
import json

import numpy as np
from numpy.typing import NDArray
from scipy.signal import find_peaks

from digitize_gridded_image import load_image

# step of the sub-pixel shift probe, in pixels; the seam lines are many
# pixels long, so the score varies smoothly on this scale
SHIFT_STEP_PX = 0.02

# a fitted seam line counts as supported when its own strength reaches this
# fraction of the median strength of all chosen seam lines of its axis
SUPPORT_FRACTION = 0.5


def seam_strengths(patch: NDArray[np.int64], axis: int) -> NDArray[np.float64]:
    """
    The colour jump between neighbouring pixels along one axis.

    The entry ``i`` holds the summed absolute RGB difference of the pixel
    pair that a seam line at position ``i + 1`` would run between. Cell
    interiors, white background, and masked areas give no jump at all.

    Parameters
    ----------
    patch
        Image patch of one panel with integer RGB channels.
    axis
        The axis along which neighbouring pixels are compared.

    Returns
    -------
    strengths
        One value per neighbour pair along the axis.
    """
    differences = np.abs(np.diff(patch, axis=axis)).astype(np.float64)
    other = tuple(i for i in range(differences.ndim) if i != axis)
    return differences.sum(axis=other)


def best_seams(
    strengths: NDArray[np.float64], cell_size: float, offset: float
) -> tuple[list[int], float, float, float]:
    """
    The integer seam positions of one panel axis, probed over one cell.

    Nominal seam positions are ``offset + cell_size * k``; a shift ``t`` in
    one full cell is probed, and each candidate's seams are rounded to
    integer pixel positions and scored by the MEAN seam strength along their
    lines.
    Averaging rather than summing keeps the score comparable between shifts
    that round a different number of lines into the window.

    Parameters
    ----------
    strengths
        Seam strength per neighbour pair, as from ``seam_strengths``.
    cell_size
        The shared cell size along this axis in pixels, its width for the
        columns and its height for the rows.
    offset
        The nominal position of the first seam inside the window.

    Returns
    -------
    seams
        The integer positions of the best seam lines, the shift that
        produced them, the score of that best lattice, and the shift of the
        strongest competing local maximum of the score curve.
    """
    length = len(strengths) + 1
    # the seam lattice continues on both sides of the nominal first seam, so
    # the integer multiples run until the seam positions leave the window
    count_left = int(np.ceil(offset / cell_size))
    count_right = int(np.ceil((length - offset) / cell_size)) + 1
    nominal = offset + cell_size * np.arange(-count_left, count_right)
    shifts = np.arange(0.0, cell_size, SHIFT_STEP_PX)
    scores = np.empty(len(shifts))
    seams_by_shift = []
    for index, shift in enumerate(shifts):
        seams = np.rint(nominal + shift).astype(int)
        seams = seams[(seams >= 1) & (seams <= length - 1)]
        seams_by_shift.append(seams)
        scores[index] = strengths[seams - 1].mean()
    best = int(np.argmax(scores))
    # the score curve is periodic in the shift, so the strongest rival phase
    # is the runner-up local maximum apart from the chosen one
    rivals, _ = find_peaks(scores)
    rival_shift = float("nan")
    for rival in rivals[np.argsort(scores[rivals])[::-1]]:
        distance = min(abs(rival - best), len(scores) - abs(rival - best))
        if distance * SHIFT_STEP_PX > 0.5:
            rival_shift = float(shifts[rival])
            break
    return (
        seams_by_shift[best].tolist(),
        float(shifts[best]),
        float(scores[best]),
        rival_shift,
    )


def supported_seams(
    seams: list[int], strengths: NDArray[np.float64]
) -> tuple[int, int]:
    """
    How many of the chosen seam lines are individually strong seams.

    A phase that a single stray strong structure can win is set apart from a
    phase the whole line family supports by comparing each chosen line's
    own strength to the median over the chosen lines.

    Parameters
    ----------
    seams
        The integer seam positions of the fitted lattice.
    strengths
        Seam strength per neighbour pair, as from ``seam_strengths``.

    Returns
    -------
    supported
        The number of chosen seam lines whose strength reaches half the
        median strength of the chosen lines, and the total number of lines.
    """
    values = strengths[np.array(seams) - 1]
    threshold = SUPPORT_FRACTION * float(np.median(values))
    supported = int(np.count_nonzero(values >= threshold))
    return supported, len(seams)


def circular_deviation(angles: NDArray[np.float64]) -> float:
    """The spread of angles in [0, 2 pi), as the max distance from the mean."""
    phases = np.exp(1j * angles)
    mean = np.angle(phases.mean())
    deviations = (angles - mean + np.pi) % (2 * np.pi) - np.pi
    return float(np.max(np.abs(deviations)))


def main() -> None:
    """Fit the seam phase of every panel and store the cell grid."""
    prs = argparse.ArgumentParser(description=__doc__ or "")
    prs.add_argument(
        "lattice_json", type=str, help="Lattice file written by panel_lattice.py."
    )
    prs.add_argument(
        "--overlay",
        action="store_true",
        help="Draw the fitted cell grid over the figure as a debug image.",
    )
    args = prs.parse_args()

    with open(args.lattice_json, encoding="utf-8") as handle:
        lattice = json.load(handle)
    img = load_image(lattice["image"]).astype(np.int64)

    cell_width = lattice["cell_size_px"]["width"]
    cell_height = lattice["cell_size_px"]["height"]
    anchor_dy = lattice["reference_box"]["y"] - lattice["outer_box"]["y0"]
    column_x0_px = lattice["column_x0_px"]
    panel_h = int(round(lattice["panel_height_px"]))
    print(f"cell sizes: {cell_width:.4f} x {cell_height:.4f} px; dy {anchor_dy}")

    results: list[dict[str, object]] = []
    phases: dict[int, list[float]] = {0: [], 1: []}
    for panel in lattice["panels"]:
        x0 = column_x0_px[panel["col"]]
        y0 = panel["y"] - anchor_dy
        # the windows of the last row and column reach past the image edge
        panel_w = lattice["column_width_px"][panel["col"]]
        patch = img[
            y0 : min(y0 + panel_h, img.shape[0]), x0 : min(x0 + panel_w, img.shape[1])
        ]
        strengths_x = seam_strengths(patch, 1)
        strengths_y = seam_strengths(patch, 0)
        # the column windows are gutter-derived and no longer hold the anchor
        # at a fixed offset, so their seam phase is nominal from the origin;
        # the row windows still start a constant anchor offset before the row
        seams_x, shift_x, score_x, rival_x = best_seams(strengths_x, cell_width, 0.0)
        seams_y, shift_y, score_y, rival_y = best_seams(
            strengths_y, cell_height, anchor_dy
        )
        supported_x, total_x = supported_seams(seams_x, strengths_x)
        supported_y, total_y = supported_seams(seams_y, strengths_y)
        # a seam's phase belongs to the global lattice at its absolute
        # position, so the window position enters when comparing panels
        phases[1].append(((x0 + seams_x[0]) % cell_width) / cell_width * 2 * np.pi)
        phases[0].append(((y0 + seams_y[0]) % cell_height) / cell_height * 2 * np.pi)
        results.append(
            {
                "row": panel["row"],
                "col": panel["col"],
                "x0": int(x0),
                "y0": int(y0),
                "width": int(panel_w),
                "seams_x": seams_x,
                "seams_y": seams_y,
                "shift_px": {"x": shift_x, "y": shift_y},
                "mean_score_px": {"x": score_x, "y": score_y},
                "rival_shift_px": {"x": rival_x, "y": rival_y},
                "supported_seams": {"x": supported_x, "y": supported_y},
                "seam_count": {"x": total_x, "y": total_y},
            }
        )
        print(
            f"panel ({panel['row']}, {panel['col']}): x {len(seams_x)} seams, "
            f"shift {shift_x:.2f} (rival {rival_x:.2f}), "
            f"{supported_x}/{total_x} lines supported | "
            f"y {len(seams_y)} seams, shift {shift_y:.2f} "
            f"(rival {rival_y:.2f}), {supported_y}/{total_y} lines supported"
        )

    # how strongly a single global seam lattice explains all panel phases
    deviation = {axis: circular_deviation(np.array(phases[axis])) for axis in (0, 1)}
    for axis, name in ((1, "x"), (0, "y")):
        cell_size = cell_width if axis == 1 else cell_height
        print(
            f"phase agreement along {name}: max deviation of the "
            f"{len(phases[axis])} panel phases from their mean is "
            f"{deviation[axis] / (2 * np.pi) * cell_size:.2f} px"
        )

    lattice["cell_grid"] = {
        "cell_size_px": {"width": cell_width, "height": cell_height},
        "phase_agreement_px": {
            "width": deviation[1] / (2 * np.pi) * cell_width,
            "height": deviation[0] / (2 * np.pi) * cell_height,
        },
        "panels": results,
    }
    with open(args.lattice_json, "w", encoding="utf-8") as handle:
        json.dump(lattice, handle, indent=2)
        handle.write("\n")
    print(f"updated {args.lattice_json}")

    if not args.overlay:
        return

    # delayed imports
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    # a seam at position s scores the pixel pair (s - 1, s), whose antialiased
    # blend pixel is row/column s - 1. Lines are drawn through the centre of
    # that seam pixel, which is the coordinate s - 0.5 in the image convention
    # of one pixel spanning [r, r + 1], so that the debug line covers the
    # visible seam instead of its border with the next cell interior.
    lines = [
        [(x0 + seam - 0.5, y0), (x0 + seam - 0.5, y0 + min(panel_h, img.shape[0] - y0))]
        for panel_grid in results
        for seam in panel_grid["seams_x"]
        for x0, y0 in [(panel_grid["x0"], panel_grid["y0"])]
    ] + [
        [
            (x0, y0 + seam - 0.5),
            (x0 + min(width, img.shape[1] - x0), y0 + seam - 0.5),
        ]
        for panel_grid in results
        for seam in panel_grid["seams_y"]
        for x0, y0, width in [(panel_grid["x0"], panel_grid["y0"], panel_grid["width"])]
    ]

    fig, view = plt.subplots(figsize=(img.shape[1] / 110, img.shape[0] / 110))
    view.imshow(
        img.astype(np.uint8),
        extent=(0, img.shape[1], img.shape[0], 0),
        interpolation="nearest",
    )
    view.add_collection(LineCollection(lines, colors="black", linewidths=0.5))
    view.set_axis_off()
    stem = args.lattice_json.rpartition(".")[0]
    filename = f"{stem}_cell_grid.png"
    fig.savefig(filename, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {filename}")


if __name__ == "__main__":
    main()
