"""Plotting wrappers for cross-context dynamic-mode comparisons.

The figures follow the manuscript TESTING notebook: mode proportions (pie charts),
reassignment Sankeys and chord diagram, representative genes per cluster pair and
violin plots of trajectory descriptors. Inputs are gene-trajectory AnnData objects
(``obs`` with ``gene``, ``cell_type``, the mode column and tradeSeq statistics).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from chronopack._core._keys import CHRONOPACK_UNS_KEY
from chronopack.pl._gene_trajectory_plots import (
    LOW_SIGNIFICANCE_COLOR,
    MODE_COLORS,
    MODE_ORDER,
    RPlotResult,
    _draw,
    _obs_table,
    _r_to_pandas,
)

LOW_SIGNIFICANCE_LABEL = "low_significance"


def _mode_table(adata, groupby: str, extra=()) -> pd.DataFrame:
    """obs columns used by the cross-context plots; the mode column is renamed as in TESTING."""
    table = _obs_table(adata, ["gene", "cell_type", groupby, "waldStat_1", *extra])
    return table.rename(columns={groupby: "broad_cluster"}) if groupby != "broad_cluster" else table


def _mode_palette(modes) -> dict[str, str]:
    palette = {m: MODE_COLORS[m] for m in modes if m in MODE_COLORS}
    if LOW_SIGNIFICANCE_LABEL in modes:
        palette[LOW_SIGNIFICANCE_LABEL] = LOW_SIGNIFICANCE_COLOR
    return palette


def mode_proportions(
    adata,
    *,
    groupby: str = "broad_cluster",
    modes=MODE_ORDER,
    colors: dict[str, str] | None = None,
    output: str | Path | None = None,
    width: float = 12,
    height: float = 5,
) -> RPlotResult:
    """Pie chart of dynamic-mode proportions, one per context (``cell_type``).

    ``data`` holds the counts and proportions per context and mode (modes in
    ``modes`` order, zero counts included; other labels are not counted).
    """
    obs = _obs_table(adata, ["cell_type", groupby])
    modes = [str(m) for m in modes]
    counts = (
        pd.crosstab(obs["cell_type"], pd.Categorical(obs[groupby].astype(str), categories=modes), dropna=False)
        .reindex(columns=modes, fill_value=0)
        .stack()
        .rename("n")
        .reset_index()
    )
    counts.columns = ["cell_type", "broad_cluster", "n"]
    counts["prop"] = counts["n"] / counts.groupby("cell_type")["n"].transform("sum")
    rendered, preview, _ = _draw(
        "chronopack_plot_mode_proportions_csv",
        output,
        tables={"data": counts},
        cluster_order=modes,
        cluster_colors=colors or {m: MODE_COLORS[m] for m in modes if m in MODE_COLORS},
        width=float(width),
        height=float(height),
    )
    return RPlotResult("mode_proportions", data=counts, output=rendered, preview_png=preview)


def mode_sankey(
    adata,
    pairs,
    *,
    groupby: str = "broad_cluster",
    top_k: int = 100,
    modes=(*MODE_ORDER, LOW_SIGNIFICANCE_LABEL),
    palette: dict[str, str] | None = None,
    treat_missing_as_low_significance: bool = True,
    low_significance_action: str = "show",
    min_flow: int = 2,
    stratum_gap: float = 3,
    stratum_width: float = 0.10,
    output: str | Path | None = None,
    width: float = 11,
    height: float = 5,
    **kwargs,
) -> RPlotResult:
    """Sankey diagrams of mode reassignment for each ``(source, target)`` context pair.

    Genes are the top ``top_k`` by ``waldStat_1`` in either context of a pair; genes
    without a mode in ``modes`` (or missing) in the other context count as
    ``low_significance``. Panels are drawn side by side (legend on the last); ``data``
    holds the flows. Extra keyword arguments go to the R
    ``plot_cluster_reassignment_sankey``.
    """
    pairs = [(str(s), str(t)) for s, t in pairs]
    table = _mode_table(adata, groupby)
    modes = [str(m) for m in modes]
    rendered, preview, value = _draw(
        "chronopack_plot_mode_sankey_csv",
        output,
        tables={"data": table},
        source_cell_types=[s for s, _ in pairs],
        target_cell_types=[t for _, t in pairs],
        top_k=int(top_k),
        treat_missing_as_low_significance=bool(treat_missing_as_low_significance),
        low_significance_action=low_significance_action,
        cluster_col="broad_cluster",
        gene_col="gene",
        rank_col="waldStat_1",
        clusters_use=modes,
        min_flow=int(min_flow),
        stratum_gap=float(stratum_gap),
        stratum_width=float(stratum_width),
        cluster_palette=palette or _mode_palette(modes),
        width=float(width),
        height=float(height),
        **kwargs,
    )
    flows = _r_to_pandas(value) if value is not None else None
    return RPlotResult(
        "mode_sankey",
        data=table if flows is None else flows,
        output=rendered,
        preview_png=preview,
        metadata={"pairs": pairs, "top_k": int(top_k)},
    )


def mode_chord(
    adata,
    contexts,
    *,
    groupby: str = "broad_cluster",
    top_k: int = 100,
    modes=(*MODE_ORDER, LOW_SIGNIFICANCE_LABEL),
    palette: dict[str, str] | None = None,
    context_colors: dict[str, str] | None = None,
    min_flow: int = 3,
    topk_set_mode: str = "union",
    treat_missing_as_low_significance: bool = True,
    show_low_significance: bool = False,
    hidden_low_significance_mode: str = "bar",
    hidden_bar_col: str = "grey75",
    hidden_bar_track_height: float = 0.04,
    celltype_track_height: float = 0.06,
    legend_inset=(-0.45, 0),
    output: str | Path | None = None,
    width: float = 8,
    height: float = 6,
    **kwargs,
) -> RPlotResult:
    """Pairwise chord diagram of mode reassignment between ``contexts``.

    Genes are the union of each context's top ``top_k`` by ``waldStat_1``; links
    of at least ``min_flow`` genes join a mode in one context to a mode in another.
    Low-significance assignments are hidden and summarised as bars by default.
    ``data`` holds the links. Extra keyword arguments go to the R
    ``plot_cluster_reassignment_chord_pairwise``.
    """
    table = _mode_table(adata, groupby)
    modes = [str(m) for m in modes]
    rendered, preview, value = _draw(
        "chronopack_plot_mode_chord_csv",
        output,
        tables={"data": table},
        compare_cell_types=[str(c) for c in contexts],
        cell_type_palette=context_colors,
        treat_missing_as_low_significance=bool(treat_missing_as_low_significance),
        celltype_track_height=float(celltype_track_height),
        cluster_col="broad_cluster",
        rank_col="waldStat_1",
        top_k=int(top_k),
        min_flow=int(min_flow),
        legend_inset=[float(x) for x in legend_inset],
        topk_set_mode=topk_set_mode,
        clusters_use=modes,
        cluster_palette=palette or _mode_palette(modes),
        show_low_significance=bool(show_low_significance),
        hidden_low_significance_mode=hidden_low_significance_mode,
        hidden_bar_col=hidden_bar_col,
        hidden_bar_track_height=float(hidden_bar_track_height),
        width=float(width),
        height=float(height),
        **kwargs,
    )
    links = _r_to_pandas(value) if value is not None else None
    return RPlotResult(
        "mode_chord",
        data=table if links is None else links,
        output=rendered,
        preview_png=preview,
        metadata={"contexts": [str(c) for c in contexts], "top_k": int(top_k)},
    )


def cluster_pair_representative_genes(
    adata,
    sce_dict: dict,
    *,
    row_context: str = "K562a",
    col_context: str = "THP1",
    clusters=("impulse-up_1", "impulse-up_2", "trend-up", "trend-down_1", "trend-down_2"),
    groupby: str = "broad_cluster",
    colours: dict[str, str] | None = None,
    palette: dict[str, str] | None = None,
    qval_cut: float = 1,
    time_map: dict[str, float] | None = None,
    offset_map: dict[str, float] | None = None,
    exclude_genes=None,
    hide_subplot_axes: bool = True,
    hide_mean_fit_points: bool = True,
    omit_points: bool = True,
    trim_before_zero: bool = True,
    top_annotation_box_rel_size: float = 0.20,
    output: str | Path | None = None,
    width: float = 10,
    height: float = 10,
    **kwargs,
) -> RPlotResult:
    """Grid of one representative gene per (row-context mode, column-context mode) pair.

    For genes in mode ``r`` in ``row_context`` and mode ``c`` in ``col_context``,
    the panel shows the highest-``waldStat_1`` gene (row context; genes in
    ``exclude_genes`` only as fallback) as fitted tradeSeq curves of both contexts.
    ``sce_dict`` maps contexts to tradeSeq SCE ``.rds`` files; ``time_map`` /
    ``offset_map`` default to the alignment stored in ``adata.uns``. ``data`` holds
    the chosen genes. Extra keyword arguments go to the R
    ``plot_cluster_pair_representative_genes`` (and its gene plot).
    """
    alignment = adata.uns.get(CHRONOPACK_UNS_KEY, {}).get("alignment", {})
    keys = [str(row_context), str(col_context)]
    time_map = time_map or alignment.get("time_map")
    offset_map = offset_map or alignment.get("offset_map")
    if time_map is None or offset_map is None:
        raise ValueError("Pass time_map and offset_map (not stored in adata.uns).")
    clusters = [str(c) for c in clusters]
    table = _mode_table(adata, groupby, extra=["qval_1"])
    rendered, preview, value = _draw(
        "chronopack_plot_cluster_pair_genes_csv",
        output,
        tables={"data": table},
        sce_paths={k: str(sce_dict[k]) for k in keys},
        cluster_list=clusters,
        row_cell_type=keys[0],
        col_cell_type=keys[1],
        cluster_col="broad_cluster",
        colours=colours or {keys[0]: "#386cb0", keys[1]: "#222021"},
        qval_cut=float(qval_cut),
        time_max_map={k: float(time_map[k]) for k in keys},
        offset_map={k: float(offset_map[k]) for k in keys},
        hide_subplot_axes=bool(hide_subplot_axes),
        hide_mean_fit_points=bool(hide_mean_fit_points),
        omit_points=bool(omit_points),
        trim_before_zero=bool(trim_before_zero),
        print_gene_table=False,
        cluster_palette=palette or _mode_palette(clusters),
        exclude_genes=None if exclude_genes is None else [str(g) for g in exclude_genes],
        top_annotation_box_rel_size=float(top_annotation_box_rel_size),
        width=float(width),
        height=float(height),
        **kwargs,
    )
    chosen = _r_to_pandas(value) if value is not None else None
    return RPlotResult(
        "cluster_pair_representative_genes",
        data=table if chosen is None else chosen,
        output=rendered,
        preview_png=preview,
        metadata={"row_context": keys[0], "col_context": keys[1]},
    )


def cluster_violin(
    adata,
    features,
    *,
    cluster_use,
    cluster_col: str = "broad_cluster",
    plot_titles=None,
    y_title: str | None = None,
    contexts=("K562a", "K562b", "THP1"),
    compare_mode: str = "across_keys_for_cluster",
    top_k: int = 100,
    topk_mode: str = "per_key",
    use_top_k: bool = True,
    trim_percentiles: bool = True,
    trim_lower: float = 0,
    trim_upper: float = 1,
    output: str | Path | None = None,
    width: float | None = None,
    height: float = 3.5,
    censored: bool = True,
    censored_color: str = "#E7298A",
    censored_label: str = "censored (at window end)",
    **kwargs,
) -> RPlotResult:
    """Violin plots of trajectory descriptors (``features``, e.g. ``t_peak``) per context.

    With the default ``compare_mode`` each panel compares ``contexts`` within mode
    ``cluster_use`` of ``cluster_col``, using each context's top ``top_k`` genes by
    ``waldStat_1``. One panel per feature, side by side (``width`` defaults to 3.5
    inches per panel); ``data`` holds the per-group summaries. Extra keyword
    arguments go to the R ``plot_cluster_violin``.

    With ``censored=True``, a feature with a ``<feature>_censored`` column in
    ``adata.obs`` (e.g. ``t_peak_censored`` from
    :func:`~chronopack.tl.compute_gene_trajectory_descriptors`) has its censored
    genes drawn as ``censored_color`` points; they stay in the violin, box and
    summaries at their stored value (the window end), and ``data`` counts them in
    ``n_censored``.
    """
    features = [features] if isinstance(features, str) else [str(f) for f in features]
    plot_titles = [""] * len(features) if plot_titles is None else [str(t) for t in plot_titles]
    if len(plot_titles) != len(features):
        raise ValueError("Pass one plot title per feature.")
    censored_cols = [
        f"{f}_censored" if censored and f"{f}_censored" in adata.obs else "" for f in features
    ]
    table = _obs_table(
        adata, ["gene", "cell_type", cluster_col, "waldStat_1", *features, *[c for c in censored_cols if c]]
    )
    rendered, preview, value = _draw(
        "chronopack_plot_cluster_violin_csv",
        output,
        tables={"data": table},
        features=features,
        plot_titles=plot_titles,
        y_title=y_title,
        print_summary=False,
        print_top_genes=False,
        keys_use=[str(c) for c in contexts],
        cluster_col=cluster_col,
        compare_mode=compare_mode,
        cluster_use=str(cluster_use),
        use_top_k=bool(use_top_k),
        top_k=int(top_k),
        topk_mode=topk_mode,
        trim_percentiles=bool(trim_percentiles),
        trim_lower=float(trim_lower),
        trim_upper=float(trim_upper),
        width=float(3.5 * len(features) if width is None else width),
        height=float(height),
        censored_cols=censored_cols,
        censored_color=censored_color,
        censored_label=censored_label,
        **kwargs,
    )
    summary = _r_to_pandas(value) if value is not None else None
    return RPlotResult(
        "cluster_violin",
        data=table if summary is None else summary,
        output=rendered,
        preview_png=preview,
        metadata={"features": features, "cluster_use": str(cluster_use)},
    )


__all__ = [
    "cluster_pair_representative_genes",
    "cluster_violin",
    "mode_chord",
    "mode_proportions",
    "mode_sankey",
]
