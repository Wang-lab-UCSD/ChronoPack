"""Matplotlib plots for cell-cell communication results."""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from chronopack.tl.c2c._helpers import to_dense


def plot_pair_overlay(
    pair_series: dict,
    pair=("TNF", "TNFRSF1A"),
    norm_mode: str = "zscore",
    nfkb_onset_time: Optional[float] = None,
    tmin_min: float = 0.0,
    tmax_min: float = 165.0,
    xtick_step: float = 30.0,
    figsize=(5, 2),
    dpi: int = 300,
    title: Optional[str] = None,
    box_aspect: Optional[float] = None,
):
    """Overlay one LR activity curve with its receiver response curve.

    ``box_aspect`` fixes the plotting area's height / width (e.g. ``390 / 445``);
    ``None`` lets it fill ``figsize``.
    """
    import matplotlib.pyplot as plt
    if pair not in pair_series:
        raise KeyError(f"{pair} not found in pair_series")

    def normalize(values):
        values = np.asarray(values, dtype=float)
        if norm_mode == "zscore":
            sd = np.nanstd(values)
            return np.zeros_like(values) if not np.isfinite(sd) or sd == 0 else (values - np.nanmean(values)) / sd
        if norm_mode == "minmax":
            rng = np.nanmax(values) - np.nanmin(values)
            return np.zeros_like(values) if not np.isfinite(rng) or rng == 0 else (values - np.nanmin(values)) / rng
        if norm_mode in {"none", "raw"}:
            return values
        raise ValueError("norm_mode must be 'zscore', 'minmax', or 'none'")

    def scale_index(index):
        arr = np.asarray(index, dtype=float)
        lo, hi = np.nanmin(arr), np.nanmax(arr)
        return np.full_like(arr, tmin_min, dtype=float) if hi == lo else (arr - lo) / (hi - lo) * (tmax_min - tmin_min) + tmin_min

    activity, response = pair_series[pair]
    x_activity = scale_index(activity.index)
    x_response = scale_index(response.index)
    y_activity = normalize(activity.to_numpy())
    y_response = normalize(response.to_numpy())

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ligand, receptor = pair
    ax.plot(x_activity, y_activity, linewidth=2.2, label=f"{ligand}->{receptor}")
    response_line, = ax.plot(x_response, y_response, linewidth=2.2, linestyle="--", label="receiver response")

    if nfkb_onset_time is not None:
        order = np.argsort(x_response)
        y_onset = np.interp(float(nfkb_onset_time), x_response[order], y_response[order])
        ax.plot(float(nfkb_onset_time), y_onset, marker="D", markersize=7, markerfacecolor="white", markeredgecolor=response_line.get_color(), markeredgewidth=1.8, linestyle="None", zorder=10)

    ax.set_xlim(tmin_min, tmax_min)
    ax.set_xticks(np.arange(tmin_min, tmax_min + 1, xtick_step))
    ax.set_xlabel("time (min)")
    ax.set_ylabel({"zscore": "z-score", "minmax": "min-max normalized"}.get(norm_mode, "value"))
    ax.set_title(title or f"{ligand}->{receptor}")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if box_aspect is not None:
        ax.set_box_aspect(box_aspect)
    plt.tight_layout()
    return fig, ax


def plot_genes_over_time(
    adata,
    gene_list: Sequence[str],
    time_key: str = "time",
    cluster_key: str = "cell_type",
    layer: Optional[str] = None,
    ncols: int = 3,
    figsize_per_panel=(4.0, 3.2),
):
    """Plot raw expression trajectories for a short list of genes."""
    import matplotlib.pyplot as plt
    time_values = pd.to_numeric(pd.Series(adata.obs[time_key]).astype(str), errors="coerce")
    cluster_values = pd.Series(adata.obs[cluster_key]).astype(str)
    mask_valid = time_values.notna()
    clusters = sorted(cluster_values[mask_valid].unique())
    colors = dict(zip(clusters, plt.rcParams["axes.prop_cycle"].by_key().get("color", [])))
    nrows = math.ceil(len(gene_list) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))
    axes = np.asarray(axes).ravel()
    matrix = adata.layers[layer] if layer is not None and layer in adata.layers else adata.X
    for ax, gene in zip(axes, gene_list):
        if gene not in adata.var_names:
            ax.text(0.5, 0.5, f"{gene}\nnot found", ha="center", va="center")
            ax.axis("off")
            continue
        expr = pd.Series(np.asarray(to_dense(matrix[:, adata.var_names.get_loc(gene)])).ravel())[mask_valid.to_numpy()].reset_index(drop=True)
        t = time_values[mask_valid].reset_index(drop=True)
        c = cluster_values[mask_valid].reset_index(drop=True)
        for cluster in clusters:
            sub = pd.DataFrame({"time": t[c == cluster], "expr": expr[c == cluster]})
            stats = sub.groupby("time").agg(mean=("expr", "mean"), n=("expr", "size"), sd=("expr", "std")).reset_index()
            stats["se"] = stats["sd"].fillna(0.0) / np.sqrt(stats["n"].clip(lower=1))
            ax.plot(stats["time"], stats["mean"], label=cluster, color=colors.get(cluster))
            ax.fill_between(stats["time"], stats["mean"] - 1.96 * stats["se"], stats["mean"] + 1.96 * stats["se"], alpha=0.2, color=colors.get(cluster))
        ax.set_title(gene)
        ax.set_xlabel("time")
        ax.set_ylabel("expression")
    for ax in axes[len(gene_list) :]:
        ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=max(1, len(labels)), frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, axes
