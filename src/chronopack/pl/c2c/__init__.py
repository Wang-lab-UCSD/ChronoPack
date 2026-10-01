"""Plotting helpers for C2C results."""

from __future__ import annotations

import numpy as np

from ._plots import plot_genes_over_time, plot_pair_overlay
from chronopack.tl.c2c import C2CResult


def interaction_scatter(
    result: C2CResult,
    *,
    p_thresh: float = 0.05,
    x_thresh: float | None = None,
    receptor_detection_clip: float | None = 0.2,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    figsize=(4.8, 4.0),
    dpi: int = 300,
    title: str = "Candidate L-R pairs",
    xlabel: str | None = None,
    ylabel: str = "Ligand temporal significance",
):
    """Scatter LR pair activity against ligand temporal significance.

    Points are colored by receptor detection in the receiver (clipped at
    ``receptor_detection_clip``). Dashed lines mark ``p_thresh`` and
    ``x_thresh``; pairs above both are labeled.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    df = result.interactions.replace([np.inf, -np.inf], np.nan).dropna(subset=["activity_sum", "logp"]).copy()
    color_col = "rec_det"
    if receptor_detection_clip is not None:
        color_col = "rec_det_clipped"
        df[color_col] = np.clip(df["rec_det"].astype(float), 0, receptor_detection_clip)
    x, y, c = df["activity_sum"], df["logp"], df[color_col]
    y_thresh = -np.log10(p_thresh)
    label_mask = (y > y_thresh) & (x >= x_thresh if x_thresh is not None else True)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.axhline(y_thresh, color="0.35", lw=2, ls="--", zorder=0)
    if x_thresh is not None:
        ax.axvline(x_thresh, color="0.35", lw=2, ls="--", zorder=0)
    sc = ax.scatter(
        x,
        y,
        s=300,
        c=c,
        cmap="Reds",
        norm=Normalize(vmin=np.nanmin(c), vmax=np.nanmax(c)),
        edgecolors="#444444",
        linewidths=np.where(y > y_thresh, 0.8, 0.6),
        alpha=0.96,
        zorder=2,
    )
    labeled = df.loc[label_mask].sort_values(["logp", "activity_sum"], ascending=[False, False])
    for i, (_, row) in enumerate(labeled.iterrows()):
        ax.annotate(
            str(row["pair"]).replace("->", r"$\rightarrow$"),
            xy=(row["activity_sum"], row["logp"]),
            xytext=(0, 7 + 3 * (i % 3)),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.78),
        )
    cbar = fig.colorbar(sc, ax=ax, pad=0.012, fraction=0.045)
    cbar.ax.set_title(color_col, fontsize=12, pad=6)
    cbar.ax.tick_params(labelsize=8)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is None:
        # Headroom above the top row of points for their labels.
        ylim = (min(float(y.min()), 0.0) - 0.35, float(y.max()) + 0.5)
    ax.set_ylim(*ylim)
    if xlabel is None:
        window = result.parameters.get("causal_window")
        xlabel = "L-R activity score" + (f" ({window[0]:g}-{window[1]:g} min)" if window else "")
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=16, pad=12)
    ax.grid(True, color="0.9", lw=0.9, alpha=0.18)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", direction="out", length=5.0, width=1.0, labelsize=8)
    ax.tick_params(axis="both", which="minor", direction="out", length=3.0, width=1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig, ax


def pair_overlay(result: C2CResult, pair: tuple[str, str], *, response=None, **kwargs):
    """Overlay one LR activity profile with a receiver response profile.

    ``response`` is a time-indexed (minutes) ``pd.Series`` of the receiver
    response to overlay. ``cp.tl.c2c.infer`` no longer computes one, so it must
    be passed unless ``result.temporal_profiles`` has a ``receiver_response``
    column.
    """
    import pandas as pd

    if result.temporal_profiles is None:
        raise ValueError("C2CResult does not contain temporal_profiles.")
    sub = result.temporal_profiles[
        (result.temporal_profiles["ligand"] == pair[0])
        & (result.temporal_profiles["receptor"] == pair[1])
    ]
    if sub.empty:
        raise KeyError(f"{pair} not found in C2CResult temporal_profiles.")
    times = sub["time"].to_numpy(float)
    activity = pd.Series(sub["lr_activity"].to_numpy(float), index=times)
    if response is None:
        if "receiver_response" not in sub.columns:
            raise ValueError("Pass response=<time-indexed receiver response pd.Series> to overlay.")
        response = pd.Series(sub["receiver_response"].to_numpy(float), index=times)
    response = pd.Series(np.asarray(response, dtype=float), index=pd.Index(response.index, dtype=float))
    common = activity.index.intersection(response.index)
    return plot_pair_overlay({pair: (activity.loc[common], response.loc[common])}, pair=pair, **kwargs)


__all__ = ["interaction_scatter", "pair_overlay", "plot_genes_over_time", "plot_pair_overlay"]
