"""Plotting helpers for C2C results."""

from __future__ import annotations

from typing import Optional

import numpy as np

from chronopack._legacy.c2c import plot_avg_vs_logp, plot_genes_over_time, plot_pair_overlay
from chronopack.tl.c2c import C2CResult


def interaction_scatter(
    result: C2CResult,
    *,
    color_col: Optional[str] = None,
    receptor_detection_clip: float | None = None,
    **kwargs,
):
    """Plot ranked LR interaction scores from a precomputed C2CResult."""
    df = result.interactions.copy()
    if receptor_detection_clip is not None:
        df["rec_det_clipped"] = np.clip(df["rec_det"], 0, receptor_detection_clip)
        color_col = color_col or "rec_det_clipped"
    return plot_avg_vs_logp(df, color_col=color_col, **kwargs)


def pair_overlay(result: C2CResult, pair: tuple[str, str], **kwargs):
    """Overlay one LR activity profile with its receiver response profile."""
    if result.temporal_profiles is None:
        raise ValueError("C2CResult does not contain temporal_profiles.")
    sub = result.temporal_profiles[
        (result.temporal_profiles["ligand"] == pair[0])
        & (result.temporal_profiles["receptor"] == pair[1])
    ]
    if sub.empty:
        raise KeyError(f"{pair} not found in C2CResult temporal_profiles.")
    import pandas as pd

    pair_series = {
        pair: (
            pd.Series(sub["lr_activity"].to_numpy(float), index=sub["time"].to_numpy(float)),
            pd.Series(sub["receiver_response"].to_numpy(float), index=sub["time"].to_numpy(float)),
        )
    }
    return plot_pair_overlay(pair_series, pair=pair, **kwargs)


__all__ = ["interaction_scatter", "pair_overlay", "plot_genes_over_time"]
