"""Composite plots of fitted tradeSeq gene-dynamics curves."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


def gene_dynamics_curves(
    curves: pd.DataFrame,
    genes,
    *,
    observations: pd.DataFrame | None = None,
    group_col: str = "lineage",
    colors: dict[str, str] | None = None,
    ncols: int = 4,
    panel_size: tuple[float, float] = (4.2, 4.0),
    box_aspect: float | None = 1.0,
    xlabel: str = "Time (min)",
    ylabel: str = "Expression",
    x_tick_step: float | None = 30,
    show_band: bool = True,
    band_alpha: float = 0.35,
    line_width: float = 3.6,
    point_size: float = 3.0,
    point_alpha: float = 0.3,
    point_edge_width: float = 0.7,
    mean_size: float = 9.4,
    mean_edge_width: float = 1.3,
    tick_label_size: float = 14,
    axis_label_size: float = 16,
    gene_title_size: float = 16.8,
    y_upper_scale: float | None = 1.4,
    y_pad: float = 0.04,
    legend_title: str | None = None,
    title: str | None = None,
    output: str | Path | None = None,
):
    """Grid of fitted NB-GAM curves, one panel per gene and one curve per group.

    ``curves`` has ``gene``, ``t``, ``mu_log1p`` (and ``lo_log1p`` / ``hi_log1p``
    for the 95% band) plus ``group_col``, e.g. stacked
    :func:`chronopack.tl.predict_gene_dynamics` outputs for the lineages of a fit or
    for separate fits. ``observations`` (stacked
    :func:`chronopack.tl.observed_gene_dynamics` outputs with the same
    ``group_col``) adds each cell's normalized expression as a point, in the
    color of its group (of its highest-weight lineage when ``assigned`` is
    given), and the per-time-point mean as a white circle: the (weighted, by
    ``weight`` when given) mean of the normalized counts, on the log1p scale.
    Panels follow ``genes`` row by row (``ncols`` per row); groups keep their order
    of appearance. Returns the matplotlib figure; ``output`` also saves it.

    The defaults reproduce the ggplot2 ``plot_tradeSeq_gene`` figures (one
    4.8 x 3.6 in plot per gene): sizes are in points as drawn (``point_size`` and
    ``mean_size`` are marker diameters including their edges; cell points have a
    ``point_edge_width`` outline in their own colour, like ggplot2's shape 19),
    ``theme_minimal`` text and grid, and axes that end at the data along time.
    ``box_aspect`` is each plotting area's height / width (1: square, x and y axes
    of equal length; ``None``: fill the panel).

    Each panel's y-axis runs from ``max(0, y_min - pad)`` to ``y_upper_scale``
    times the panel's highest fitted mean, where ``y_min`` is the lowest observed
    value, per-time-point mean or band lower bound (fitted mean without a band)
    and ``pad`` is ``y_pad`` of the span from ``y_min`` to the upper limit. Cells
    above the upper limit are cut off. ``y_upper_scale=None`` keeps matplotlib's
    automatic limits.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, MultipleLocator

    genes = [str(g) for g in genes]
    missing = sorted(set(genes) - set(curves["gene"].astype(str)))
    if missing:
        raise KeyError(f"No fitted curves for: {', '.join(missing)}")
    groups = list(dict.fromkeys(curves[group_col].astype(str)))
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = {g: (colors or {}).get(g, palette[i % len(palette)]) for i, g in enumerate(groups)}

    means = None
    if observations is not None:
        obs = observations.assign(gene=observations["gene"].astype(str), **{group_col: observations[group_col].astype(str)})
        weight = obs["weight"] if "weight" in obs else pd.Series(1.0, index=obs.index)
        obs = obs.assign(_w=weight, _wy=weight * np.expm1(obs["expr_log1p"]))
        sums = obs.groupby(["gene", group_col, "t"], sort=True)[["_w", "_wy"]].sum().reset_index()
        means = sums[sums["_w"] > 0].assign(mean_log1p=lambda d: np.log1p(d["_wy"] / d["_w"]))
        if "assigned" in obs:
            obs = obs[obs["assigned"].astype(bool)]

    # Marker areas (pt^2) from drawn diameters; a marker edge is centred on its outline.
    point_area = max(point_size - point_edge_width, 0.0) ** 2
    mean_area = max(mean_size - mean_edge_width, 0.0) ** 2

    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_size[0] * ncols, panel_size[1] * nrows), squeeze=False)
    for ax, gene in zip(axes.flat, genes):
        sub = curves[curves["gene"].astype(str) == gene]
        t_values = [sub["t"]]
        if observations is not None:
            pts = obs[obs["gene"] == gene]
            t_values.append(pts["t"])
            for group in groups:
                g = pts[pts[group_col] == group]
                ax.scatter(g["t"], g["expr_log1p"], s=point_area, color=colors[group], alpha=point_alpha,
                           lw=point_edge_width, rasterized=True, zorder=1)
        for group in groups:
            g = sub[sub[group_col].astype(str) == group].sort_values("t")
            if g.empty:
                continue
            if show_band and {"lo_log1p", "hi_log1p"}.issubset(g.columns):
                ax.fill_between(g["t"], g["lo_log1p"], g["hi_log1p"], color=colors[group], alpha=band_alpha,
                                lw=0, zorder=2)
            ax.plot(g["t"], g["mu_log1p"], color=colors[group], lw=line_width, label=group, zorder=3,
                    solid_capstyle="butt")
        if means is not None:
            for group in groups:
                m = means[(means["gene"] == gene) & (means[group_col] == group)]
                ax.scatter(m["t"], m["mean_log1p"], s=mean_area, facecolor="white", edgecolor=colors[group],
                           lw=mean_edge_width, zorder=4)
        ax.set_title(gene, fontsize=gene_title_size, fontweight="bold", loc="left")
        if box_aspect is not None:
            ax.set_box_aspect(box_aspect)
        if x_tick_step:
            ax.xaxis.set_major_locator(MultipleLocator(x_tick_step))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))  # like scales::pretty_breaks
        t_all = pd.concat(t_values).to_numpy(float)
        ax.set_xlim(np.nanmin(t_all), np.nanmax(t_all))  # like coord_cartesian(expand = FALSE)
        if y_upper_scale is not None:
            groups_sub = sub[sub[group_col].astype(str).isin(groups)]
            lows = [groups_sub["lo_log1p" if show_band and "lo_log1p" in groups_sub else "mu_log1p"]]
            if observations is not None:
                lows += [pts["expr_log1p"], means.loc[means["gene"] == gene, "mean_log1p"]]
            y_min = float(np.nanmin(pd.concat(lows).to_numpy(float)))
            upper = y_upper_scale * float(np.nanmax(groups_sub["mu_log1p"].to_numpy(float)))
            if np.isfinite(upper) and upper > y_min:
                ax.set_ylim(max(0.0, y_min - y_pad * (upper - y_min)), upper)
        ax.grid(True, which="major", color="#e0e0e0", lw=1.4, zorder=0)  # theme_minimal: grey88
        ax.set_axisbelow(True)
        ax.spines[:].set_visible(False)
        ax.tick_params(labelsize=tick_label_size, length=0, pad=5, colors="black")
    for ax in axes.flat[len(genes):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel(xlabel, fontsize=axis_label_size)
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel, fontsize=axis_label_size)
    if title:
        fig.suptitle(title, fontsize=gene_title_size, fontweight="bold")
    fig.tight_layout()
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
               title=legend_title, fontsize=tick_label_size * 0.8, title_fontsize=tick_label_size)
    if output is not None:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight", dpi=300)
    return fig


__all__ = ["gene_dynamics_curves"]
