"""
Measure the shared cell width and height of a multi-panel figure.

Each panel plots its data as small patches (cells) of a constant pixel size,
separated by thin antialiased seams, so the seam period shows up as a peak in
the Fourier spectrum of the panel's brightness differences along each axis.
A single panel's spectrum is dominated by the large-scale structure of the
disk it shows, so the amplitude spectra of all panels are stacked: the cell
period is common to every panel and survives, while the disk structure
averages partly away. This script reports the stacked peak period for the two
axes, with a sub-pixel estimate from parabolic peak interpolation and the
spread of the per-panel peaks, and plots the stacked spectra as the input for
fitting the integer cell grid inside the panels.
"""

import argparse
import json

import numpy as np
from numpy.typing import NDArray

from digitize_gridded_image import load_image

# plausible cell period in pixels; the cell seams in these figures sit near
# 5 px, while the structure of the disk dominates periods above ~50 px
MIN_PERIOD_PX = 3.0
MAX_PERIOD_PX = 9.0


def difference_signal(panel: NDArray[np.int64], axis: int) -> NDArray[np.float64]:
    """
    The mean absolute brightness difference along one axis of a panel.

    Averaging the luminance first and the absolute differences over the other
    axis afterwards keeps each row (or column) of seams as an independent
    sample of the seam pattern.

    Parameters
    ----------
    panel
        Image patch of one panel with integer RGB channels.
    axis
        The axis along which neighbouring pixels are differenced (1 for the
        horizontal seam profile averaged over rows, 0 for the vertical one
        averaged over columns).

    Returns
    -------
    signal
        One value per neighbour pair along the axis.
    """
    luminance = panel.mean(axis=2)
    differences = np.abs(np.diff(luminance, axis=axis))
    return differences.mean(axis=1 - axis)


def refine_peak(values: NDArray[np.float64], index: int) -> float:
    """The sub-bin peak position by parabolic interpolation."""
    left, middle, right = values[index - 1], values[index], values[index + 1]
    return index + 0.5 * (left - right) / (left - 2 * middle + right)


def peak_period(spectrum: NDArray[np.float64], length: int) -> float:
    """
    The period of the strongest peak of a spectrum inside the cell band.

    Parameters
    ----------
    spectrum
        Amplitude spectrum of a difference signal of length ``length``.
    length
        The length of the differenced signal the spectrum was computed from.

    Returns
    -------
    period
        The peak period in pixels, refined to sub-bin accuracy.
    """
    frequencies = np.fft.rfftfreq(length)
    band = (frequencies > 1 / MAX_PERIOD_PX) & (frequencies < 1 / MIN_PERIOD_PX)
    index = int(np.argmax(np.where(band, spectrum, 0.0)))
    return length / refine_peak(spectrum, index)


def main() -> None:
    """Stack the panel spectra and report the cell period candidates."""
    prs = argparse.ArgumentParser(description=__doc__ or "")
    prs.add_argument(
        "lattice_json", type=str, help="Lattice file written by panel_lattice.py."
    )
    args = prs.parse_args()

    with open(args.lattice_json, encoding="utf-8") as handle:
        lattice = json.load(handle)
    img = load_image(lattice["image"]).astype(np.int64)

    # a panel window starts at its column's gutter-derived origin and at the
    # measured row anchor shifted back by the anchor's offset in the window
    anchor_dy = lattice["reference_box"]["y"] - lattice["outer_box"]["y0"]
    column_x0_px = lattice["column_x0_px"]
    panel_h = int(round(lattice["panel_height_px"]))

    titles = {1: "horizontal seams (cell width)", 0: "vertical seams (cell height)"}
    signals: dict[int, list[NDArray[np.float64]]] = {0: [], 1: []}
    for panel in lattice["panels"]:
        x0 = column_x0_px[panel["col"]]
        y0 = panel["y"] - anchor_dy
        # the windows of the last row and column reach past the image edge
        # and are clipped, so only the shortest per-axis signal is stacked
        panel_w = lattice["column_width_px"][panel["col"]]
        patch = img[
            y0 : min(y0 + panel_h, img.shape[0]), x0 : min(x0 + panel_w, img.shape[1])
        ]
        for axis in (0, 1):
            signals[axis].append(difference_signal(patch, axis))

    # crop every signal to the shortest one of its axis so that the stacked
    # spectra share one frequency grid
    spectrums: dict[int, list[NDArray[np.float64]]] = {0: [], 1: []}
    lengths: dict[int, int] = {
        axis: min(len(signal) for signal in signals[axis]) for axis in (0, 1)
    }
    periods_px: dict[int, list[float]] = {0: [], 1: []}
    for axis in (0, 1):
        for signal in signals[axis]:
            cropped = signal[: lengths[axis]]
            spectrum = np.abs(np.fft.rfft(cropped - cropped.mean()))
            spectrums[axis].append(spectrum)
            periods_px[axis].append(peak_period(spectrum, lengths[axis]))

    # delayed imports
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stem = args.lattice_json.rpartition(".")[0]
    cell_size_px: dict[str, float] = {}
    figure, views = plt.subplots(1, 2, figsize=(11, 4))
    for axis, view in zip((1, 0), views):
        # the signals share the common cropped length, so their spectra share
        # one frequency grid and the amplitudes can be averaged bin by bin
        stacked = np.mean(np.stack(spectrums[axis]), axis=0)
        length = lengths[axis]
        frequencies = np.fft.rfftfreq(length)
        periods = np.divide(
            1.0,
            frequencies,
            out=np.zeros_like(frequencies),
            where=frequencies > 0,
        )
        band = (periods > MIN_PERIOD_PX) & (periods < MAX_PERIOD_PX)
        view.plot(periods[band], stacked[band])
        view.set_xlabel("period (px)")
        view.set_ylabel("stacked amplitude")
        view.set_title(f"{titles[axis]}\n{stem}", fontsize=10)
        cell_size = peak_period(stacked, length)
        cell_size_px["width" if axis == 1 else "height"] = cell_size
        print(
            f"{titles[axis]}: stacked peak {cell_size:.4f} px "
            f"(per-panel peaks: mean {np.mean(periods_px[axis]):.4f}, "
            f"min {np.min(periods_px[axis]):.4f}, "
            f"max {np.max(periods_px[axis]):.4f})"
        )
    figure.tight_layout()
    filename = f"{stem}_cell_fft.png"
    figure.savefig(filename, dpi=110)
    plt.close(figure)
    print(f"wrote {filename}")

    # store the measured cell sizes in the lattice so the seam-phase fit
    # chains from them instead of repeating the estimate
    lattice["cell_size_px"] = {
        "width": cell_size_px["width"],
        "height": cell_size_px["height"],
    }
    with open(args.lattice_json, "w", encoding="utf-8") as handle:
        json.dump(lattice, handle, indent=2)
        handle.write("\n")
    print(f"updated {args.lattice_json}")


if __name__ == "__main__":
    main()
