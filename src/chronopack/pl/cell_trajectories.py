"""Plotting helpers for cell-state trajectory and fate outputs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from chronopack._core._keys import (
    CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    OT_GRAPH_KEY,
    OT_UMAP_KEY,
    TERMINAL_FATE_PROBS_KEY,
)
from chronopack.tl._trajectories import (
    PlotHardClusterMassProportionsOverTime,
    fate_umap_2x2_scanpy,
    mass_retention_plots_and_umap,
)


def fate_probabilities(
    adata,
    *,
    fate_key="terminal_fate_probs_traj",
    fate_names_key="terminal_fate_probs_traj_fates",
    embedding_basis=OT_UMAP_KEY,
    fate_obs_prefix="fate",
    which_two=(0, 1),
    max_entropy_norm=True,
    size=15,
):
    """Plot two fate probabilities, fate entropy and max fate probability.

    Also writes per-fate probabilities plus ``{fate_obs_prefix}_entropy`` and
    ``{fate_obs_prefix}_maxprob`` to ``adata.obs``, as the workflow has always
    done; these columns are saved with the trajectory outputs.
    """
    return fate_umap_2x2_scanpy(
        adata,
        fate_key=fate_key,
        fate_names_key=fate_names_key,
        umap_basis=embedding_basis,
        fate_obs_prefix=fate_obs_prefix,
        which_two=which_two,
        max_entropy_norm=max_entropy_norm,
        size=size,
    )


def proportions_over_time(adata, *, cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY, time_key="time"):
    """Plot hard cell-state trajectory-cluster proportions over time."""
    return PlotHardClusterMassProportionsOverTime(
        adata,
        cluster_key=cluster_key,
        time_key=time_key,
    )


def mass_retention(
    adata,
    *,
    S=None,
    cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    time_key="time",
    embedding_basis=OT_UMAP_KEY,
    **kwargs,
):
    """Plot cell-state trajectory mass retention from a precomputed square graph."""
    if S is None:
        if OT_GRAPH_KEY not in adata.obsp:
            raise KeyError(f"Missing adata.obsp[{OT_GRAPH_KEY!r}]; pass S explicitly or run cp.tl.get_trajectories().")
        S = adata.obsp[OT_GRAPH_KEY]
    return mass_retention_plots_and_umap(
        adata,
        S,
        cluster_key=cluster_key,
        time_key=time_key,
        umap_basis=embedding_basis,
        **kwargs,
    )


def _save(fig, output, dpi):
    if output is not None:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")


def _style_embedding_axis(ax, legend_title=""):
    ax.set_xlabel("UMAP1", fontsize=14)
    ax.set_ylabel("UMAP2", fontsize=14)
    ax.set_title("")
    leg = ax.get_legend()
    if leg is not None:
        leg.set_title(legend_title, prop={"size": 13})
        for txt in leg.get_texts():
            txt.set_fontsize(12)


def time_minutes(
    adata,
    *,
    embedding_basis=OT_UMAP_KEY,
    time_key="time",
    minutes_per_tag=10,
    first_tag=1,
    obs_key="Time (min)",
    size=40,
    figsize=(4, 4),
    output=None,
    dpi=300,
):
    """Embedding coloured by sampling time in minutes.

    Time tags become ``(tag - first_tag) * minutes_per_tag`` minutes, stored in
    ``adata.obs[obs_key]`` as an ordered categorical with the latest time first
    (top of the legend). Returns the figure; ``output`` also saves it.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    t_idx = pd.Series(adata.obs[time_key], index=adata.obs_names).astype(str).astype(float).astype(int)
    t_min = (t_idx - int(first_tag)) * minutes_per_tag
    cats = np.sort(t_min.unique())
    adata.obs[obs_key] = pd.Categorical(t_min, categories=cats, ordered=True)
    adata.obs[obs_key] = adata.obs[obs_key].cat.reorder_categories(list(cats[::-1]), ordered=True)

    with plt.rc_context({"figure.figsize": figsize}):
        axs = sc.pl.embedding(adata, basis=embedding_basis, color=[obs_key], show=False, size=size)
    axs = np.ravel([axs]) if not isinstance(axs, (list, np.ndarray)) else np.ravel(axs)
    for ax in axs:
        _style_embedding_axis(ax, legend_title=obs_key)
    fig = axs[0].figure
    _save(fig, output, dpi)
    return fig


def trajectory_cell_types(
    adata,
    *,
    cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    labels=None,
    palette=("#2C6DA4", "#B84A4A"),
    obs_key="traj_cell_type",
    embedding_basis=OT_UMAP_KEY,
    size=40,
    figsize=(4, 4),
    output=None,
    dpi=300,
):
    """Embedding coloured by named trajectory clusters.

    ``labels`` maps cluster ids to names (legend order follows the mapping;
    default: the sorted ids). The names are stored in ``adata.obs[obs_key]``
    with ``palette`` as their colours. Returns the figure.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    clusters = adata.obs[cluster_key].astype(str)
    if labels is None:
        labels = {c: c for c in sorted(clusters.unique())}
    labels = {str(k): v for k, v in labels.items()}
    order = list(dict.fromkeys(labels.values()))
    if len(palette) < len(order):
        raise ValueError(f"palette has {len(palette)} colours for {len(order)} labels.")
    adata.obs[obs_key] = pd.Categorical(clusters.map(labels), categories=order, ordered=True)
    palette = list(palette)[: len(order)]
    adata.uns[f"{obs_key}_colors"] = palette

    with plt.rc_context({"figure.figsize": figsize}):
        axs = sc.pl.embedding(adata, basis=embedding_basis, color=obs_key, palette=palette, show=False, size=size)
    ax = axs if not isinstance(axs, (list, np.ndarray)) else np.ravel(axs)[0]
    _style_embedding_axis(ax, legend_title="")
    fig = ax.figure
    _save(fig, output, dpi)
    return fig


def fate_probability_panels(
    adata,
    *,
    fate_key=TERMINAL_FATE_PROBS_KEY,
    fate_names_key=f"{TERMINAL_FATE_PROBS_KEY}_fates",
    embedding_basis=OT_UMAP_KEY,
    fate_obs_prefix="fate",
    which_two=(0, 1),
    max_entropy_norm=True,
    size=15,
    figsize=(5, 4.5),
    output=None,
    dpi=300,
    save_panels=False,
    outdir=None,
    file_prefix="fate_umap",
    panel_figsize=(5, 4.5),
    rasterize_points=True,
    cbar_ticklabel_size=12,
    cbar_tick_width=1.2,
    cbar_tick_length=4.0,
    umap_border_width=1.0,
):
    """Publication-style 2x2 of two fate probabilities, fate entropy and max fate probability.

    Same panels and ``adata.obs`` columns as :func:`fate_probabilities`
    (``fate[which_two[0]]`` in Reds, ``fate[which_two[1]]`` in Blues), without
    titles, axis labels or ticks, and with styled colour bars. ``output`` saves
    the 2x2 figure; ``save_panels=True`` also saves each panel as PDF and PNG
    (points rasterized) to ``outdir``. Returns the figure.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    def _safe(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")

    assert fate_key in adata.obsm, f"Need adata.obsm['{fate_key}']"
    assert fate_names_key in adata.uns, f"Need adata.uns['{fate_names_key}']"

    F = np.asarray(adata.obsm[fate_key])
    fates = list(adata.uns[fate_names_key])

    assert F.shape[0] == adata.n_obs
    assert F.shape[1] == len(fates)

    i0, i1 = which_two
    if not (0 <= i0 < len(fates) and 0 <= i1 < len(fates)):
        raise ValueError(f"which_two={which_two} out of range for {len(fates)} fates")

    for j, name in enumerate(fates):
        k = f"{fate_obs_prefix}_{name}"
        if k not in adata.obs:
            adata.obs[k] = F[:, j].astype(np.float32)

    eps = 1e-12
    P = np.clip(F, eps, 1.0)
    ent = -np.sum(P * np.log(P), axis=1)
    if max_entropy_norm and len(fates) > 1:
        ent = ent / np.log(len(fates))
    adata.obs[f"{fate_obs_prefix}_entropy"] = ent.astype(np.float32)
    adata.obs[f"{fate_obs_prefix}_maxprob"] = np.max(F, axis=1).astype(np.float32)

    k0 = f"{fate_obs_prefix}_{fates[i0]}"
    k1 = f"{fate_obs_prefix}_{fates[i1]}"
    kH = f"{fate_obs_prefix}_entropy"
    kM = f"{fate_obs_prefix}_maxprob"

    def _style_umap_axis(ax):
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(umap_border_width)

    def _style_colorbar_axis(cax):
        cax.set_ylabel("")
        cax.tick_params(labelsize=cbar_ticklabel_size, width=cbar_tick_width, length=cbar_tick_length)
        for spine in cax.spines.values():
            spine.set_linewidth(cbar_tick_width)

    def _clean_axes_and_colorbars(fig, main_axes):
        main_axes = list(main_axes)
        for ax in main_axes:
            _style_umap_axis(ax)
        # Scanpy colorbars are extra axes in the figure
        for a in fig.axes:
            if a not in main_axes:
                _style_colorbar_axis(a)

    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)

    sc.pl.embedding(adata, basis=embedding_basis, color=k0, cmap="Reds", vmin=0, vmax=1,
                    size=size, ax=axes[0, 0], show=False, title="")
    sc.pl.embedding(adata, basis=embedding_basis, color=k1, cmap="Blues", vmin=0, vmax=1,
                    size=size, ax=axes[0, 1], show=False, title="")
    sc.pl.embedding(adata, basis=embedding_basis, color=kH, cmap="viridis", vmin=0,
                    vmax=1 if max_entropy_norm else None,
                    size=size, ax=axes[1, 0], show=False, title="")
    sc.pl.embedding(adata, basis=embedding_basis, color=kM, cmap="viridis", vmin=0, vmax=1,
                    size=size, ax=axes[1, 1], show=False, title="")

    _clean_axes_and_colorbars(fig, axes.ravel())
    _save(fig, output, dpi)

    if save_panels:
        if outdir is None:
            raise ValueError("save_panels=True needs outdir.")
        Path(outdir).mkdir(parents=True, exist_ok=True)
        panels = [
            (k0, "Reds", 0, 1, f"{file_prefix}__prob__{_safe(fates[i0])}"),
            (k1, "Blues", 0, 1, f"{file_prefix}__prob__{_safe(fates[i1])}"),
            (kH, "viridis", 0, (1 if max_entropy_norm else None), f"{file_prefix}__entropy"),
            (kM, "viridis", 0, 1, f"{file_prefix}__maxprob"),
        ]
        for color_key, cmap, vmin, vmax, base_name in panels:
            fig1, ax1 = plt.subplots(1, 1, figsize=panel_figsize, constrained_layout=True)
            sc.pl.embedding(adata, basis=embedding_basis, color=color_key, cmap=cmap,
                            vmin=vmin, vmax=vmax, size=size, ax=ax1, show=False, title="")
            if rasterize_points:
                for coll in ax1.collections:
                    coll.set_rasterized(True)
            _clean_axes_and_colorbars(fig1, [ax1])
            for ext in (".pdf", ".png"):
                fig1.savefig(Path(outdir) / (base_name + ext), dpi=dpi, bbox_inches="tight")
            plt.close(fig1)

    return fig


def trajectory_endpoints(
    adata,
    *,
    cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    labels=None,
    palette=None,
    endpoint_key="cell_type",
    time_key="time",
    endpoint_times=None,
    minutes_per_tag=10,
    first_tag=1,
    endpoint_darken=1.0,
    label_offset=0.015,
    other_color="#B0B0B0",
    embedding_basis=OT_UMAP_KEY,
    point_size=15,
    figsize=(5, 5),
    legend_ncol=None,
    output=None,
    dpi=300,
):
    """Embedding coloured by named trajectory clusters, with the first and last time labelled.

    THP-1 + K562 CreateFigures endpoint UMAP. ``labels`` maps cluster ids to names
    (legend order, above the plot, follows the mapping; default: the sorted ids) and
    ``palette`` maps names (or lists, in label order) to colours; unlisted clusters
    are drawn in ``other_color``. All points have the same size, cells of the
    ``endpoint_times`` (default: the first and last time tag) on top, coloured
    ``endpoint_darken`` times their group colour (1 keeps it). For each value of
    ``obs[endpoint_key]`` the median position of its cells at each endpoint is
    labelled ``(tag - first_tag) * minutes_per_tag`` min, shifted down (first) or
    up (last) by ``label_offset`` of the y range. Returns the figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex, to_rgb
    from matplotlib.lines import Line2D

    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    if labels is None:
        labels = {c: c for c in sorted(set(clusters))}
    labels = {str(k): v for k, v in labels.items()}
    order = list(dict.fromkeys(labels.values()))
    if palette is None:
        palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    if not isinstance(palette, dict):
        if len(palette) < len(order):
            raise ValueError(f"palette has {len(palette)} colours for {len(order)} labels.")
        palette = dict(zip(order, palette))
    missing = [name for name in order if name not in palette]
    if missing:
        raise KeyError(f"No colour for label(s): {', '.join(missing)}")

    t = pd.to_numeric(pd.Series(adata.obs[time_key]).astype(str)).to_numpy()
    if endpoint_times is None:
        endpoint_times = (np.nanmin(t), np.nanmax(t))
    endpoint_times = list(endpoint_times)
    groups = np.array([labels.get(c) for c in clusters], dtype=object)

    def _darken(color):
        return to_hex(np.clip(np.array(to_rgb(color)) * endpoint_darken, 0, 1))

    is_end = np.isin(t, endpoint_times)
    colors = np.array([
        other_color if g is None else (_darken(palette[g]) if end else palette[g])
        for g, end in zip(groups, is_end)
    ])

    coords = np.asarray(adata.obsm[embedding_basis])[:, :2]
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(coords[~is_end, 0], coords[~is_end, 1], c=colors[~is_end], s=point_size,
               linewidths=0, alpha=0.85, zorder=1)
    ax.scatter(coords[is_end, 0], coords[is_end, 1], c=colors[is_end], s=point_size,
               linewidths=0, alpha=0.95, zorder=2)
    ax.set_xlabel("UMAP1", fontsize=14)
    ax.set_ylabel("UMAP2", fontsize=14)
    ax.set_title("")
    ax.set_xticks([])
    ax.set_yticks([])

    y_range = coords[:, 1].max() - coords[:, 1].min()
    endpoint_groups = adata.obs[endpoint_key].astype(str).to_numpy()
    first, last = min(endpoint_times), max(endpoint_times)
    for value in dict.fromkeys(endpoint_groups):
        for tag in endpoint_times:
            pts = coords[(endpoint_groups == value) & (t == tag)]
            if pts.shape[0] == 0:
                continue
            dy = -label_offset if tag == first else (label_offset if tag == last else 0.0)
            ax.text(
                np.median(pts[:, 0]), np.median(pts[:, 1]) + dy * y_range,
                f"{(tag - first_tag) * minutes_per_tag:g} min",
                fontsize=11, fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="black", alpha=0.9), zorder=5,
            )

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[name], markersize=9, label=name)
        for name in order
    ]
    ax.legend(handles=handles, frameon=False, fontsize=12, loc="lower center", bbox_to_anchor=(0.5, 1.02),
              ncol=len(order) if legend_ncol is None else legend_ncol, columnspacing=0.2,
              handletextpad=0.2, borderaxespad=0.0)
    fig.subplots_adjust(top=0.86)
    _save(fig, output, dpi)
    return fig


# --- two-gene fate scatter / UMAP helpers (CreateFigures.ipynb) ---

def _volcano_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("volcano", ["black", "darkred", "red", "orange", "yellow", "white"])


def _coerce_cmap(cmap):
    """``"volcano"`` (black-red-yellow-white) as a colormap; other names pass through."""
    if isinstance(cmap, str) and cmap.lower() == "volcano":
        return _volcano_cmap()
    return cmap


def _resolve_two_gene_param(param, gene_left, gene_right, name="param"):
    """Scalar, length-2 sequence (left, right) or dict keyed by gene."""
    if isinstance(param, dict):
        return (
            _coerce_cmap(param.get(gene_left, "viridis")),
            _coerce_cmap(param.get(gene_right, "viridis")),
        )
    if isinstance(param, (list, tuple)):
        if len(param) != 2:
            raise ValueError(f"{name} must have length 2 if provided as a list/tuple.")
        return _coerce_cmap(param[0]), _coerce_cmap(param[1])
    return _coerce_cmap(param), _coerce_cmap(param)


def _first_time_bin(adata, time_key="time"):
    t = adata.obs[time_key]
    try:
        tn = pd.to_numeric(t, errors="coerce")
        if np.isfinite(tn).any():
            return np.nanmin(tn)
    except Exception:
        pass
    return sorted(t.unique(), key=lambda x: str(x))[0]


def _dense_column(v):
    return v.toarray().ravel() if sp.issparse(v) else np.asarray(v).ravel()


def _pull_expr(adata, gene, mask, layer_preference=("log1p",)):
    if adata.raw is not None and gene in adata.raw.var_names:
        return _dense_column(adata.raw[:, [gene]].X)[mask]
    for lyr in layer_preference:
        if lyr in getattr(adata, "layers", {}) and gene in adata.var_names:
            return _dense_column(adata.layers[lyr][:, adata.var_names.get_loc(gene)])[mask]
    if gene in adata.var_names:
        return _dense_column(adata.X[:, adata.var_names.get_loc(gene)])[mask]
    return np.full(int(mask.sum()), np.nan)


def _gene_use_raw(adata, gene):
    if adata.raw is not None and gene in adata.raw.var_names:
        return True
    if gene in adata.var_names:
        return False
    raise ValueError(f"{gene} not found in adata.var_names or adata.raw.var_names")


def _style_fate_scatter_ax(
    ax, gene, x, y, color, dot_size=110, x_label="", y_label="Log Expression",
    font_title=28, font_tick=22, font_annot=24, xlim=(0, 1), xticks=(0, 0.5, 1.0),
):
    from scipy.stats import spearmanr

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)

    ax.scatter(x[valid], y[valid], s=dot_size, color=color, alpha=0.95, edgecolors="white", linewidths=0.5)

    if valid.sum() >= 2 and np.unique(x[valid]).size > 1 and np.unique(y[valid]).size > 1:
        rho, _ = spearmanr(x[valid], y[valid])
    else:
        rho = np.nan
    rho_label = f"r = {rho:.2f}" if np.isfinite(rho) else "r = NA"
    ax.text(
        0.98, 0.98, rho_label, transform=ax.transAxes, ha="right", va="top",
        fontsize=font_annot, color="black", zorder=10,
        bbox=dict(facecolor="white", edgecolor="black", linewidth=1.2, alpha=0.95, boxstyle="round,pad=0.30"),
    )

    ax.set_title(gene, fontsize=font_title, fontweight="bold", pad=8)
    ax.set_xlabel(x_label, fontsize=font_title)
    ax.set_ylabel(y_label, fontsize=font_title)
    ax.set_xlim(*xlim)
    ax.set_xticks(list(xticks))
    ax.set_xticklabels([f"{value:g}" for value in xticks])
    ax.grid(color="0.92", linewidth=1)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.6)
    ax.spines["bottom"].set_linewidth(1.6)
    ax.tick_params(axis="both", labelsize=font_tick, length=5, width=1.4, direction="out")
    ax.margins(x=0.00, y=0.12)
    return rho


def _new_axes_after_plot(fig, axes_before):
    """Axes a Scanpy call added to the figure (its colorbars)."""
    return [ax for ax in fig.axes if ax not in axes_before]


def _style_two_gene_umap_axis(ax, font_title=28, font_tick=22):
    """No title; UMAP2 label on the right."""
    ax.set_title("")
    ax.set_xlabel("UMAP1", fontsize=font_title)
    ax.set_ylabel("UMAP2", fontsize=font_title, rotation=270, labelpad=28)
    ax.yaxis.set_label_position("right")
    ax.tick_params(axis="both", labelsize=font_tick, length=5, width=1.4, direction="out")


def _place_colorbar_left(
    cax, umap_ax, label="Log Expression", tick_interval=1.0, tick_fontsize=18,
    label_fontsize=22, width=0.012, pad=0.018,
):
    """Move a vertical colorbar to the left of its UMAP axis (width/pad in figure coordinates)."""
    from matplotlib.ticker import FormatStrFormatter, MultipleLocator

    umap_pos = umap_ax.get_position()
    cax.set_position([umap_pos.x0 - pad - width, umap_pos.y0, width, umap_pos.height])
    cax.yaxis.set_major_locator(MultipleLocator(tick_interval))
    cax.yaxis.set_major_formatter(FormatStrFormatter("%d"))
    cax.yaxis.set_ticks_position("left")
    cax.yaxis.set_label_position("left")
    cax.tick_params(axis="y", labelsize=tick_fontsize, length=4, width=1.2, direction="out",
                    labelleft=True, labelright=False)
    cax.set_ylabel(label, fontsize=label_fontsize, rotation=90, labelpad=18)


def two_gene_fate_scatter(
    adata,
    genes,
    *,
    colors=("#CB4856", "#4F6FAE"),
    time_key="time",
    early_bin=None,
    fate_key="fate_0",
    embedding_basis=OT_UMAP_KEY,
    layer_preference=("log1p",),
    figsize=(11, 10),
    dot_size=110,
    umap_dot_size=40,
    cmap="viridis",
    expression_vmax=None,
    x_label=None,
    y_label="Log Expression",
    xlim=(0, 1),
    xticks=(0, 0.5, 1.0),
    font_title=28,
    font_tick=22,
    font_annot=24,
    row_hspace=0.35,
    col_wspace=0.35,
    figure_left=0.12,
    figure_right=0.97,
    figure_bottom=0.08,
    figure_top=0.97,
    colorbar_label="Log Expression",
    colorbar_tick_interval=1.0,
    colorbar_tick_font_scale=1.0,
    colorbar_label_font_scale=1.0,
    colorbar_width=0.012,
    colorbar_pad=0.018,
    output=None,
    dpi=300,
):
    """Two genes against a fate probability in the earliest time bin, plus their UMAPs.

    Top row: expression (``.raw``, else the first of ``layer_preference``, else
    ``.X``) of each gene versus ``obs[fate_key]`` for cells in ``early_bin``
    (default: the smallest ``obs[time_key]``), with the Spearman r. Bottom
    row: each gene on ``embedding_basis`` (``.X``, vmin 0, vmax per gene from
    ``expression_vmax``, default ``"p99"``), colour bars on the left.
    ``cmap`` may be one name (``"volcano"`` included), a pair or a dict by gene.
    The Spearman r values are in ``fig.spearman_r``. Returns the figure.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    if len(genes) != 2:
        raise ValueError("genes must contain exactly two gene names.")
    if len(colors) != 2:
        raise ValueError("colors must contain exactly two colors.")

    gene_left, gene_right = genes
    color_left, color_right = colors
    cmap_left, cmap_right = _resolve_two_gene_param(cmap, gene_left, gene_right, name="cmap")
    x_label = fate_key if x_label is None else x_label

    if expression_vmax is None:
        expression_vmax = {}
    if isinstance(expression_vmax, dict):
        vmax_left = expression_vmax.get(gene_left, "p99")
        vmax_right = expression_vmax.get(gene_right, "p99")
    else:
        vmax_left = vmax_right = expression_vmax

    if early_bin is None:
        early_bin = _first_time_bin(adata, time_key)
    mask = adata.obs[time_key].isin([early_bin]).to_numpy()
    if not np.any(mask):
        raise ValueError(f"No cells in early_bin={early_bin!r} for obs[{time_key!r}]")
    if fate_key not in adata.obs.columns:
        raise KeyError(f"{fate_key!r} is not present in adata.obs.")
    x = pd.to_numeric(adata.obs.loc[mask, fate_key], errors="coerce").to_numpy()

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    layout = dict(left=figure_left, right=figure_right, bottom=figure_bottom, top=figure_top,
                  hspace=row_hspace, wspace=col_wspace)
    fig.subplots_adjust(**layout)

    scatter_style = dict(dot_size=dot_size, x_label=x_label, y_label=y_label, font_title=font_title,
                         font_tick=font_tick, font_annot=font_annot, xlim=xlim, xticks=xticks)
    rho = {}
    colorbars = []
    for col, (gene, color) in enumerate([(gene_left, color_left), (gene_right, color_right)]):
        y = _pull_expr(adata, gene, mask, layer_preference=layer_preference)
        rho[gene] = _style_fate_scatter_ax(axes[0, col], gene, x, y, color, **scatter_style)

    for col, (gene, gene_cmap, vmax) in enumerate([(gene_left, cmap_left, vmax_left), (gene_right, cmap_right, vmax_right)]):
        axes_before = set(fig.axes)
        sc.pl.embedding(adata, basis=embedding_basis, color=gene, use_raw=_gene_use_raw(adata, gene),
                        cmap=gene_cmap, size=umap_dot_size, vmin=0, vmax=vmax, ax=axes[1, col], show=False)
        colorbars.append(_new_axes_after_plot(fig, axes_before))
        _style_two_gene_umap_axis(axes[1, col], font_title=font_title, font_tick=font_tick)

    # Scanpy may alter the layout
    fig.subplots_adjust(**layout)
    for col, cbar_axes in enumerate(colorbars):
        if len(cbar_axes) > 0:
            _place_colorbar_left(
                cbar_axes[0], axes[1, col], label=colorbar_label, tick_interval=colorbar_tick_interval,
                tick_fontsize=font_tick * colorbar_tick_font_scale,
                label_fontsize=font_title * colorbar_label_font_scale,
                width=colorbar_width, pad=colorbar_pad,
            )

    fig.spearman_r = rho
    _save(fig, output, dpi)
    return fig


def top_fate_correlations(
    res,
    top_n=20,
    q_threshold=0.05,
    *,
    neg_dot="#2C6DA4",
    pos_dot="#B84A4A",
    label_fontsize=8,
    point_size=12,
    figsize=(8, 8.0),
    x_step=2.0,
    gap_steps=3.0,
    left_pad_steps=0.8,
    right_pad_steps=0.75,
    label_offset_frac=0.06,
    panel_margin_frac=0.05,
    y_pad_frac=0.16,
    tick_label_size=10,
    axis_label_size=12,
    abs_neg_ticks=False,
    top_ylabel="Fate1 r",
    bot_ylabel="Fate0 r",
    xlabel="Gene",
    show_corner_breaks=False,
    output=None,
    dpi=300,
):
    """Top positively (upper panel) and negatively (lower panel) fate-correlated genes.

    ``res`` has ``gene``, ``r`` and ``qval`` columns, e.g. from
    :func:`chronopack.tl.early_fate_correlations`; genes with
    ``qval < q_threshold`` are ranked by r. ``abs_neg_ticks`` shows |r| on the
    lower panel. The plotted genes are in ``fig.top_genes`` (``pos``, ``neg``).
    Returns the figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    sig = res.loc[(res["qval"] < q_threshold) & (~res["r"].isna()), ["gene", "r"]].copy()

    neg = sig.sort_values("r", ascending=True).head(top_n).reset_index(drop=True)
    pos = sig.sort_values("r", ascending=False).head(top_n).reset_index(drop=True)

    neg["x"] = np.arange(len(neg)) * x_step
    start_pos = ((len(neg) - 1) + gap_steps) * x_step if len(neg) else 0.0
    pos["x"] = np.arange(len(pos)) * x_step + start_pos

    def _lims(d):
        if d.empty:
            return (-1, -0.8)
        lo, hi = d["r"].min(), d["r"].max()
        pad = max(y_pad_frac, label_offset_frac + panel_margin_frac + 0.04) * max(1e-9, hi - lo)
        return (lo - pad, hi + pad)

    ylim_pos = _lims(pos)
    ylim_neg = _lims(neg)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, figsize=figsize, gridspec_kw=dict(height_ratios=[1, 1], hspace=0.05)
    )
    ax_top.set_ylim(*ylim_pos)
    ax_bot.set_ylim(*ylim_neg)

    def draw(ax, df, dot_color, ylim):
        if df.empty:
            return
        ax.scatter(df["x"], df["r"], s=point_size, color=dot_color, zorder=3)
        yrng = ylim[1] - ylim[0]
        yoff = label_offset_frac * yrng
        ycap_top = ylim[1] - panel_margin_frac * yrng
        for x, y, g in zip(df["x"], df["r"], df["gene"]):
            ylab = min(y + yoff, ycap_top)
            ax.text(x, ylab, g, rotation=90, ha="center", va="bottom",
                    fontsize=label_fontsize, color="black", clip_on=False, zorder=4)

    draw(ax_top, pos, pos_dot, ylim_pos)
    draw(ax_bot, neg, neg_dot, ylim_neg)

    for ax in (ax_top, ax_bot):
        ax.grid(axis="y", color="0.93", lw=1)
        ax.set_xticks([])
        ax.tick_params(labelsize=tick_label_size, length=4, direction="out")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.spines["left"].set_linewidth(1.6)
        ax.spines["bottom"].set_linewidth(1.6)

    if abs_neg_ticks:
        ax_bot.yaxis.set_major_formatter(FuncFormatter(lambda v, pos=None: f"{abs(v):.1f}"))

    xmax = max(neg["x"].max() if len(neg) else 0, pos["x"].max() if len(pos) else 0)
    ax_top.set_xlim(-left_pad_steps * x_step, xmax + right_pad_steps * x_step)

    ax_top.set_ylabel(top_ylabel, fontsize=axis_label_size)
    ax_bot.set_ylabel(bot_ylabel, fontsize=axis_label_size)
    ax_bot.set_xlabel(xlabel, fontsize=axis_label_size)

    if show_corner_breaks:
        d = 0.012
        kw_t = dict(transform=ax_top.transAxes, color="k", clip_on=False, lw=1.2)
        kw_b = dict(transform=ax_bot.transAxes, color="k", clip_on=False, lw=1.2)
        ax_top.plot((-d, +d), (1 - d, 1 + d), **kw_t)
        ax_top.plot((1 - d, 1 + d), (1 - d, 1 + d), **kw_t)
        ax_bot.plot((-d, +d), (-d, +d), **kw_b)
        ax_bot.plot((1 - d, 1 + d), (-d, +d), **kw_b)

    fig.subplots_adjust(top=0.98, bottom=0.12)
    fig.top_genes = {"pos": pos["gene"].tolist(), "neg": neg["gene"].tolist()}
    _save(fig, output, dpi)
    return fig


__all__ = [
    "fate_probabilities",
    "fate_probability_panels",
    "mass_retention",
    "proportions_over_time",
    "time_minutes",
    "top_fate_correlations",
    "trajectory_cell_types",
    "trajectory_endpoints",
    "two_gene_fate_scatter",
]
