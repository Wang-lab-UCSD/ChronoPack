"""Small subset of :mod:`skimage.filters` needed by Scanpy Scrublet.

The maintained ChronoPack environment used for verification does not include
scikit-image. Scanpy only needs ``threshold_minimum`` when Scrublet is allowed
to choose doublet thresholds automatically, so this shim avoids mutating the
conda environment while preserving that code path.
"""

from __future__ import annotations

import numpy as np


def _local_maxima(values: np.ndarray) -> np.ndarray:
    return np.flatnonzero((values[1:-1] > values[:-2]) & (values[1:-1] > values[2:])) + 1


def _uniform_filter1d_width3(values: np.ndarray) -> np.ndarray:
    padded = np.pad(values, 1, mode="edge")
    return (padded[:-2] + padded[1:-1] + padded[2:]) / 3.0


def threshold_minimum(image=None, nbins: int = 256, max_num_iter: int = 10000, *, hist=None):
    """Return the minimum between two histogram peaks.

    This mirrors the scikit-image algorithm used by Scrublet closely enough for
    the one-dimensional simulated-doublet score distributions used here.
    """
    if hist is not None:
        counts, bin_centers = hist
        counts = np.asarray(counts, dtype=float)
        bin_centers = np.asarray(bin_centers, dtype=float)
    else:
        if image is None:
            raise ValueError("Either image or hist must be provided.")
        arr = np.asarray(image, dtype=float).ravel()
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            raise ValueError("Cannot threshold an empty image.")
        counts, edges = np.histogram(arr, bins=nbins)
        counts = counts.astype(float, copy=False)
        bin_centers = (edges[:-1] + edges[1:]) / 2.0

    smooth = counts.astype(float, copy=True)
    maxima = _local_maxima(smooth)
    for _ in range(max_num_iter):
        smooth = _uniform_filter1d_width3(smooth)
        maxima = _local_maxima(smooth)
        if len(maxima) < 3:
            break

    if len(maxima) != 2:
        raise RuntimeError("Unable to find two maxima in histogram.")

    lo, hi = maxima
    threshold_idx = lo + int(np.argmin(smooth[lo : hi + 1]))
    return float(bin_centers[threshold_idx])
