"""Plotting wrappers for gene-expression trajectory AnnData outputs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core._keys import GENE_TRAJECTORY_CLUSTERING_KEY, GENE_TRAJECTORY_KEY


@dataclass
class RPlotResult:
    """Small result object for ChronoPack plots."""

    plot_type: str
    data: pd.DataFrame
    output: str | None = None
    backend: str = "R/ggplot2"
    metadata: dict[str, Any] | None = None
    figure: Any | None = None
    preview_png: bytes | None = field(default=None, repr=False)

    def __repr__(self) -> str:
        return f"RPlotResult(plot_type={self.plot_type!r}, output={self.output!r})"

    # Inline display width relative to the rendered image (high-dpi PNGs look too large at 1:1).
    display_scale: float = field(default=0.5, repr=False)

    # Jupyter shows the plot inline from a PNG preview: unlike SVG, PNG outputs are
    # displayed even in notebooks Jupyter does not trust (e.g. executed by a script).
    def _repr_png_(self):
        png = self.preview_png
        if png is None and self.output and Path(self.output).suffix.lower() == ".png" and Path(self.output).exists():
            png = Path(self.output).read_bytes()
        if png is None:
            return None
        width_px = int.from_bytes(png[16:20], "big")  # PNG IHDR width
        return png, {"width": max(1, int(width_px * self.display_scale))}


def _trajectory_matrix(adata) -> np.ndarray:
    return np.asarray(adata.X, dtype=float)


def _time_grid(adata, n_timepoints: int) -> np.ndarray:
    if "time" in adata.var:
        out = np.asarray(adata.var["time"], dtype=float)
    else:
        out = np.arange(n_timepoints, dtype=float)
    if out.shape[0] != n_timepoints:
        raise ValueError(f"Stored time grid has length {out.shape[0]}, expected {n_timepoints}.")
    return out


def _gene_labels(adata) -> pd.DataFrame:
    labels = pd.DataFrame(index=adata.obs_names.astype(str))
    labels["row_id"] = labels.index
    labels["gene"] = (
        adata.obs["gene"].astype(str).to_numpy()
        if "gene" in adata.obs
        else labels["row_id"].str.replace(r"_(?=[^_]+$).*", "", regex=True)
    )
    labels["cell_type"] = adata.obs["cell_type"].astype(str).to_numpy() if "cell_type" in adata.obs else ""
    return labels


def _select_rows(adata, *, genes=None) -> pd.Index:
    obs_names = pd.Index(adata.obs_names.astype(str))
    if genes is None:
        return obs_names
    requested = pd.Index([str(g) for g in genes])
    direct = requested[requested.isin(obs_names)]
    selected = list(direct)
    remaining = requested.difference(direct)
    if len(remaining):
        labels = _gene_labels(adata)
        selected.extend(labels.index[labels["gene"].isin(remaining)].tolist())
    selected_idx = pd.Index(selected).drop_duplicates()
    found_genes = pd.Index(_gene_labels(adata).loc[selected_idx, "gene"]) if len(selected_idx) else pd.Index([])
    missing = requested.difference(selected_idx.union(found_genes))
    if len(selected_idx) == 0 or len(missing) > 0:
        raise KeyError(f"Could not find requested genes/rows: {', '.join(missing.astype(str))}")
    return selected_idx


def _z_row(matrix: np.ndarray) -> np.ndarray:
    mu = np.nanmean(matrix, axis=1, keepdims=True)
    sd = np.nanstd(matrix, axis=1, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return (matrix - mu) / sd


def _trajectory_long(adata, *, rows: pd.Index, normalize: bool) -> pd.DataFrame:
    fit = _trajectory_matrix(adata)
    row_lookup = pd.Series(np.arange(adata.n_obs), index=adata.obs_names.astype(str))
    idx = row_lookup.loc[rows].to_numpy()
    values = fit[idx, :]
    if normalize:
        values = _z_row(values)
    time = _time_grid(adata, values.shape[1])

    meta = _gene_labels(adata).loc[rows].reset_index(drop=True)
    for col in ["waldStat_1", "peak_log2fc", "t_onset", "t_peak", "broad_cluster", "cluster_fine_induced", "cluster_fine_repressed"]:
        if col in adata.obs:
            meta[col] = adata.obs.loc[rows, col].to_numpy()

    out = pd.DataFrame(values, index=rows, columns=[f"t{i}" for i in range(values.shape[1])])
    out.insert(0, "row_id", rows.astype(str))
    out = out.melt(id_vars="row_id", var_name="time_index", value_name="value")
    out["time_index"] = out["time_index"].str.replace("t", "", regex=False).astype(int)
    out["time"] = time[out["time_index"].to_numpy()]
    return out.merge(meta, on="row_id", how="left")


# Default palettes of TESTING.ipynb: `cluster_colors` (named dynamic modes, also the R
# chronopack_mode_palette), the grey of `low_significance`, and `cat_colors` (numbered clusters).
MODE_ORDER = ("impulse-up_1", "impulse-up_2", "trend-up", "transient-down", "trend-down_1", "trend-down_2")
MODE_COLORS = {
    "impulse-up_1": "#F8766D",
    "impulse-up_2": "#D84A42",
    "trend-up": "#C29A33",
    "transient-down": "#8B78C6",
    "trend-down_1": "#619CFF",
    "trend-down_2": "#3C6FB6",
}
LOW_SIGNIFICANCE_COLOR = "#D1D1D1"
CATEGORICAL_COLORS = (
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#8C564B", "#6A3D9A", "#66A61E",
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999", "#1B9E77", "#7570B3",
    "#E7298A", "#A6761D", "#D95F02", "#17BECF", "#BCBD22", "#AEC7E8", "#FF9896", "#98DF8A", "#C5B0D5", "#C49C94",
)


def _default_colors(labels) -> dict[str, str] | list[str]:
    """``MODE_COLORS`` when every label is a named dynamic mode, else ``CATEGORICAL_COLORS`` in label order."""
    levels = pd.unique(pd.Series(labels).dropna().astype(str))
    if set(levels) <= set(MODE_COLORS):
        return dict(MODE_COLORS)
    n = max(60, len(levels))
    return [CATEGORICAL_COLORS[i % len(CATEGORICAL_COLORS)] for i in range(n)]


def _to_r(value):
    """Convert a Python argument to an R vector (dicts become named vectors)."""
    from rpy2.robjects import BoolVector, FloatVector, IntVector, StrVector

    if isinstance(value, dict):
        vec = _to_r(list(value.values()))
        vec.names = StrVector([str(k) for k in value])
        return vec
    if isinstance(value, (bool, np.bool_)):
        return BoolVector([bool(value)])
    if isinstance(value, (int, np.integer)):
        return IntVector([int(value)])
    if isinstance(value, (float, np.floating)):
        return FloatVector([float(value)])
    if isinstance(value, (str, Path)):
        return StrVector([str(value)])
    if isinstance(value, (list, tuple, np.ndarray, pd.Index, pd.Series)):
        values = list(value)
        if all(isinstance(v, (bool, np.bool_)) for v in values) and values:
            return BoolVector([bool(v) for v in values])
        if all(isinstance(v, (int, float, np.integer, np.floating)) for v in values) and values:
            return FloatVector([float(v) for v in values])
        return StrVector([str(v) for v in values])
    return value


def _r_to_pandas(value) -> pd.DataFrame | None:
    from rpy2.robjects import conversion, default_converter, pandas2ri
    from rpy2.robjects.vectors import DataFrame as RDataFrame

    if not isinstance(value, RDataFrame):
        return None
    with conversion.localconverter(default_converter + pandas2ri.converter):
        out = conversion.rpy2py(value)
    return out.reset_index(drop=True)


def _draw(
    plot_function: str,
    output,
    *,
    tables: dict[str, pd.DataFrame],
    resource: str = "dimred_plotting_helpers.R",
    **kwargs,
):
    """Draw R ``plot_function`` to ``output`` and a PNG preview in one call.

    Each table is passed as a CSV file (argument ``<name>_path``); other arguments are
    converted with ``_to_r`` (``None`` leaves the R default). Returns the output path,
    the preview PNG bytes and the R return value; nothing is drawn without ``output``.
    """
    if output is None:
        return None, None, None
    output = str(Path(output))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:
        raise ImportError("R-backed plotting requires optional dependency 'rpy2'.") from exc
    source_r_resource(resource)
    with TemporaryDirectory() as tmpdir:
        args = {}
        for name, table in tables.items():
            path = Path(tmpdir) / f"{name}.csv"
            table.to_csv(path, index=False)
            args[f"{name}_path"] = str(path)
        outputs = [output]
        if Path(output).suffix.lower() != ".png":
            outputs.append(str(Path(tmpdir) / "preview.png"))
        args.update({k: _to_r(v) for k, v in kwargs.items() if v is not None})
        value = r[plot_function](output=StrVector(outputs), **args)
        preview = Path(outputs[-1]).read_bytes() if Path(outputs[-1]).exists() else None
    return output, preview, value


def _render(plot_function: str, data: pd.DataFrame, output, **kwargs) -> tuple[str | None, bytes | None]:
    """Draw ``plot_function`` (R) to ``output``; return the path and a PNG preview for inline display."""
    rendered, preview, _ = _draw(plot_function, output, tables={"data": data}, **kwargs)
    return rendered, preview


def gene_trajectory(
    adata,
    genes,
    *,
    fit_key: str = GENE_TRAJECTORY_KEY,
    clustering_key: str | None = GENE_TRAJECTORY_CLUSTERING_KEY,
    normalize: bool = False,
    output: str | Path | None = None,
    width: float = 9,
    height: float = 6,
    title: str | None = None,
) -> RPlotResult:
    """Plot fitted trajectories for selected genes or gene-context rows."""
    del fit_key, clustering_key
    rows = _select_rows(adata, genes=genes)
    data = _trajectory_long(adata, rows=rows, normalize=normalize)
    rendered, preview = _render(
        "chronopack_plot_gene_trajectory_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
        normalize=bool(normalize),
    )
    return RPlotResult("gene_trajectory", data=data, output=rendered, preview_png=preview, metadata={"normalize": normalize})


def _obs_table(adata, columns) -> pd.DataFrame:
    """``row_id`` plus the requested ``adata.obs`` columns, as plain values."""
    missing = [c for c in columns if c not in adata.obs]
    if missing:
        raise KeyError(f"Trajectory observation column(s) not available: {', '.join(missing)}.")
    out = pd.DataFrame({"row_id": adata.obs_names.astype(str)})
    for col in dict.fromkeys(columns):
        values = adata.obs[col]
        out[col] = values.astype(str).to_numpy() if isinstance(values.dtype, pd.CategoricalDtype) else values.to_numpy()
    return out


def trajectory_features(
    adata,
    *,
    color: str = "broad_cluster",
    basis: str | None = None,
    x: str = "PC1",
    y: str = "PC2",
    cat_colors=None,
    output: str | Path | None = None,
    width: float = 6.5,
    height: float = 5,
    title: str | None = None,
    censored: bool = True,
    censored_color: str = "#E7298A",
    censored_label: str = "censored (at window end)",
    **kwargs,
) -> RPlotResult:
    """Scatter of trajectories in a feature space, coloured by ``color`` (TESTING ``feature_plot_pc_simple``).

    ``x`` and ``y`` are ``traj.obs`` columns (descriptors), or with ``basis`` the
    components (``"PC1"``, ``"PC2"``, ...) of ``traj.obsm[basis]``, e.g. the ``X_pca``
    coordinates stored by the clustering notebook; trajectories without coordinates
    there (NaN) are left out. Categorical colours default to the dynamic-mode palette
    (``MODE_COLORS``) when every label is a named mode. Extra keyword arguments go to
    the R ``feature_plot_pc_simple``.

    With ``censored=True``, a numeric ``color`` with a ``<color>_censored`` column
    (e.g. ``t_peak_censored``) draws its censored trajectories in ``censored_color``,
    labelled ``censored_label`` in the legend, instead of the colour scale's top
    colour; the scale covers the uncensored values only.
    """
    data = _gene_labels(adata).reset_index(drop=True)
    censored_col = f"{color}_censored"
    use_censored = (
        censored
        and censored_col in adata.obs
        and pd.api.types.is_numeric_dtype(adata.obs[color])
        and censored_col not in (x, y)
    )
    if use_censored:
        data[censored_col] = adata.obs[censored_col].astype("boolean").fillna(False).astype(bool).to_numpy()
    for name in dict.fromkeys([x, y, color] if basis is None else [color]):
        if name not in adata.obs:
            raise KeyError(f"Trajectory observation column '{name}' is not available.")
        values = adata.obs[name]
        data[name] = values.astype(str).to_numpy() if isinstance(values.dtype, pd.CategoricalDtype) else values.to_numpy()
    if basis is not None:
        if basis not in adata.obsm:
            raise KeyError(f"Trajectory embedding adata.obsm['{basis}'] is not available.")
        coords = np.asarray(adata.obsm[basis], dtype=float)
        for name in (x, y):
            data[name] = coords[:, int(str(name).removeprefix("PC")) - 1]
        data = data.dropna(subset=[x, y]).reset_index(drop=True)
    categorical = not pd.api.types.is_numeric_dtype(data[color])
    if cat_colors is None and categorical:
        cat_colors = _default_colors(data[color])
    rendered, preview, _ = _draw(
        "chronopack_plot_trajectory_features_csv",
        output,
        tables={"data": data},
        feature=str(color),
        x=str(x),
        y=str(y),
        cat_colors=cat_colors,
        title=title,
        width=float(width),
        height=float(height),
        **({"censored_col": censored_col, "censored_color": censored_color, "censored_label": censored_label}
           if use_censored else {}),
        **kwargs,
    )
    return RPlotResult("trajectory_features", data=data, output=rendered, preview_png=preview, metadata={"x": x, "y": y, "color": color})


def cluster_medoids_and_top_genes(
    adata,
    groupby="broad_cluster",
    *,
    features="X_cluster_features",
    score_key: str = "waldStat_1",
    top_n: int = 9,
    gene_selection: str = "wald",
    seed: int = 13,
    show_median_iqr: bool = False,
    show_sil_text: bool = True,
    use_silhouette_alpha: bool = True,
    show_gene_titles: bool = False,
    gene_title_mode: str = "rownames",
    cat_colors=None,
    cluster_order=None,
    title: str | None = None,
    output: str | Path | None = None,
    width: float = 13,
    height: float | None = None,
    **kwargs,
) -> RPlotResult:
    """One row per cluster: medoid trajectory and top genes (TESTING ``plot_cluster_medoids_and_top_genes``).

    Panels show row z-scored trajectories (``adata.X``) over ``adata.var["time"]``;
    with ``show_median_iqr`` the first panel is the cluster median and IQR. Medoids
    and silhouettes use Euclidean distances between the clustering ``features``: an
    ``adata.obsm`` key (default, as stored by the clustering notebook), a DataFrame
    indexed by ``obs_names`` (e.g. ``GeneTrajectoryClusteringResult.features``) or an
    array aligned with ``adata``. ``groupby`` is an ``obs`` column or labels aligned
    with ``adata`` (Series indexed by ``obs_names``, e.g. candidate assignments).
    Genes are ranked by ``obs[score_key]``. ``title=None`` keeps the R default title;
    ``""`` drops it. ``height`` defaults to the number of clusters + 1 (inches).
    Extra keyword arguments go to the R function.
    """
    obs_names = pd.Index(adata.obs_names.astype(str))
    if isinstance(groupby, str):
        if groupby not in adata.obs:
            raise KeyError(f"Trajectory observation column '{groupby}' is not available.")
        labels = pd.Series(adata.obs[groupby].to_numpy(), index=obs_names)
    elif isinstance(groupby, pd.Series):
        labels = groupby.copy()
        labels.index = labels.index.astype(str)
        labels = labels.reindex(obs_names)
    else:
        labels = pd.Series(np.asarray(groupby), index=obs_names)
    labels = labels.astype(object).where(labels.notna(), None)
    if labels.isna().any():
        raise ValueError("Every trajectory needs a cluster label; subset adata first.")
    numeric = pd.to_numeric(labels, errors="coerce")
    if numeric.notna().all() and (numeric % 1 == 0).all():
        labels = numeric.astype(int)  # numbered clusters sort numerically in R
    else:
        labels = labels.astype(str)

    if isinstance(features, str):
        if features not in adata.obsm:
            raise KeyError(f"Clustering features adata.obsm['{features}'] are not available.")
        feat = pd.DataFrame(np.asarray(adata.obsm[features], dtype=float), index=obs_names)
    elif isinstance(features, pd.DataFrame):
        feat = features.copy()
        feat.index = feat.index.astype(str)
        feat = feat.reindex(obs_names)
    else:
        feat = pd.DataFrame(np.asarray(features, dtype=float), index=obs_names)
    if feat.isna().to_numpy().any():
        raise ValueError("Clustering features are missing for some trajectories.")
    feat.columns = [f"f{i}" for i in range(feat.shape[1])]

    if score_key not in adata.obs:
        raise KeyError(f"Trajectory observation column '{score_key}' is not available.")
    z = pd.DataFrame(_z_row(_trajectory_matrix(adata)), index=obs_names, columns=adata.var_names.astype(str))
    obs = pd.DataFrame({
        "row_id": obs_names,
        "cluster": labels.to_numpy(),
        "score": pd.to_numeric(adata.obs[score_key], errors="coerce").to_numpy(),
    })
    n_clusters = labels.nunique()
    if cat_colors is None:
        cat_colors = _default_colors(labels)
    rendered, preview, _ = _draw(
        "chronopack_plot_cluster_medoids_csv",
        output,
        tables={
            "traj": z.reset_index(names="row_id"),
            "features": feat.reset_index(names="row_id"),
            "obs": obs,
            "time": pd.DataFrame({"time": _time_grid(adata, z.shape[1])}),
        },
        top_n=int(top_n),
        gene_selection=gene_selection,
        seed=int(seed),
        show_median_iqr=bool(show_median_iqr),
        show_sil_text=bool(show_sil_text),
        use_silhouette_alpha=bool(use_silhouette_alpha),
        show_gene_titles=bool(show_gene_titles),
        gene_title_mode=gene_title_mode,
        cat_colors=cat_colors,
        cluster_order=None if cluster_order is None else [str(c) for c in cluster_order],
        title=title,
        width=float(width),
        height=float(n_clusters + 1 if height is None else height),
        **kwargs,
    )
    return RPlotResult(
        "cluster_medoids_and_top_genes",
        data=obs,
        output=rendered,
        preview_png=preview,
        metadata={"n_clusters": int(n_clusters), "top_n": int(top_n)},
    )


def clustering_metrics(
    clustering_result,
    *,
    output: str | Path | None = None,
    metrics=("WSS", "p90_centroid_distance", "resolution", "mean_silhouette"),
    show=("k_elbow", "k_silhouette"),
    width: float = 9,
    height: float = 6,
) -> RPlotResult:
    """Plot candidate-k clustering diagnostics without selecting k."""
    import matplotlib.pyplot as plt

    data = clustering_result.metrics.copy()
    if "k" not in data and "K" in data:
        data["k"] = data["K"]
    fig, axes = plt.subplots(len(metrics), 1, figsize=(width, height), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, metric in zip(axes, metrics):
        ax.plot(data["k"], data[metric], marker="o", linewidth=1.2, markersize=3)
        if "k_elbow" in show and metric == "elbow_score" or ("k_elbow" in show and metric in {"WSS", "resolution"}):
            k_elbow = int(data.loc[data["elbow_score"].idxmax(), "k"])
            ax.axvline(k_elbow, color="#D55E00", linestyle="--", linewidth=1)
        if "k_silhouette" in show and metric == "mean_silhouette":
            k_sil = int(data.loc[data["mean_silhouette"].idxmax(), "k"])
            ax.scatter([k_sil], [data.loc[data["mean_silhouette"].idxmax(), metric]], color="#0072B2", s=42, zorder=3)
        ax.set_ylabel(metric)
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("k")
    fig.tight_layout()
    rendered = None
    if output is not None:
        output = str(Path(output))
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150, bbox_inches="tight")
        rendered = output
    return RPlotResult("clustering_metrics", data=data, output=rendered, backend="matplotlib", figure=fig)


VOLCANO_COLUMNS = ("signedMaxLog2FC_start", "qval_1", "waldStat_1")


def fdr_wald_break(asso_df, fdr_thresh: float = 0.01) -> float:
    """``waldStat_1`` of the gene whose ``qval_1`` is closest to ``fdr_thresh`` (first on ties)."""
    q = pd.to_numeric(asso_df["qval_1"], errors="coerce").to_numpy(float)
    w = pd.to_numeric(asso_df["waldStat_1"], errors="coerce").to_numpy(float)
    ok = np.isfinite(q) & np.isfinite(w)
    if not ok.any():
        raise ValueError("No genes with finite qval_1 and waldStat_1.")
    return float(w[ok][np.argmin(np.abs(q[ok] - float(fdr_thresh)))])


def association_volcano(
    asso_df: pd.DataFrame,
    *,
    fdr_thresh: float = 0.01,
    l2fc_thresh: float = 0.5,
    top_labels: int = 12,
    y_break: float | str | None = "fdr",
    output: str | Path | None = None,
    width: float = 7.5,
    height: float = 6.5,
    dpi: int = 500,
    **kwargs,
) -> RPlotResult:
    """Volcano of a single-lineage tradeSeq association table (TradeSeqFits ``plot_association_volcano``).

    ``asso_df`` is indexed by gene with ``waldStat_1``, ``qval_1`` (e.g.
    ``asso_df.1trajectory.nknots4.csv``) and ``signedMaxLog2FC_start``
    (:func:`chronopack.tl.signed_max_log2fc_from_start`). The y-axis is log10
    Wald, compressed below ``y_break`` (Wald units) by ``y_lower_scale``:
    ``"fdr"`` (default) uses the Wald of the gene whose ``qval_1`` is closest to
    ``fdr_thresh`` (also where the dashed line is drawn), a number is used as is,
    and ``None`` uses the ``y_break_quantile`` of the non-significant Walds.
    Points with ``qval_1 < fdr_thresh`` are coloured by sign when
    ``|log2FC| >= l2fc_thresh``; the ``top_labels`` genes by Wald are labelled.
    Labels are stacked just above their points (``label_x_nudge`` /
    ``label_y_offset_frac`` / ``label_min_vgap_frac``), then overlapping ones are
    nudged apart with ggrepel (``stack_resolve_overlaps=False`` turns that off);
    ``label_method="repel"`` lays all labels out with ggrepel instead.
    Extra keyword arguments (``x_lim``, ``y_lim``, ``y_lower_scale``,
    ``label_*``, ``axis_*_size``, ...) go to the R function. ``metadata`` holds
    the ``y_break`` used and, once drawn, the labelled genes.
    """
    missing = [c for c in VOLCANO_COLUMNS if c not in asso_df]
    if missing:
        raise KeyError(f"Association table is missing column(s): {', '.join(missing)}.")
    data = asso_df.loc[:, list(VOLCANO_COLUMNS)].apply(pd.to_numeric, errors="coerce")
    data.insert(0, "row_id", asso_df.index.astype(str))
    data = data.reset_index(drop=True)
    if isinstance(y_break, str):
        if y_break != "fdr":
            raise ValueError("y_break must be a number, 'fdr' or None.")
        y_break = fdr_wald_break(data, fdr_thresh)
    rendered, preview, value = _draw(
        "chronopack_plot_association_volcano_csv",
        output,
        tables={"data": data},
        fdr_thresh=float(fdr_thresh),
        l2fc_thresh=float(l2fc_thresh),
        top_labels=int(top_labels),
        y_break=None if y_break is None else float(y_break),
        width=float(width),
        height=float(height),
        dpi=int(dpi),
        **kwargs,
    )
    labels = _r_to_pandas(value) if value is not None else None
    return RPlotResult(
        "association_volcano",
        data=data,
        output=rendered,
        preview_png=preview,
        metadata={"y_break": y_break, "fdr_thresh": float(fdr_thresh), "n_genes": int(len(data)), "labels": labels},
    )



RESPONSE_SPECIFICITY_COLUMNS = ("waldStat_1", "waldStat_2", "qval", "log2FC_t1", "log2FC_t2")


def response_specificity_volcano(
    combined: pd.DataFrame,
    conditions=("Condition A", "Condition B"),
    *,
    output: str | Path | None = None,
    width: float = 10,
    height: float = 8.75,
    dpi: int = 300,
    **kwargs,
) -> RPlotResult:
    """Two-context response-specificity volcano (VolcanoPlot ``plot_association_volcano``).

    ``combined`` is indexed by gene with ``waldStat_1`` / ``waldStat_2``, ``qval`` and
    ``log2FC_t1`` / ``log2FC_t2`` (context 1 = ``conditions[0]``), e.g. from
    :func:`chronopack.tl.combine_context_associations`. The x-axis is
    ``|log2FC_t1| - |log2FC_t2|``, the y-axis log10 of the larger Wald statistic
    (``wald_mode="difference"``: their difference), compressed below ``y_break``
    (Wald units; default: the ``y_break_quantile`` of the Walds with
    ``qval >= fdr_thresh``) by ``y_lower_scale``. A gene is specific to a context when
    its |log2FC| there is >= ``dyn_thresh_condA`` / ``dyn_thresh_condB`` and <=
    ``flat_thresh_*`` in the other (and ``qval < fdr_thresh`` if ``require_fdr``).
    Labels: ``label_top_wald`` genes by Wald plus the ``label_top_condA`` /
    ``label_top_condB`` most specific genes (else ``top_labels``). Extra keyword
    arguments (thresholds, ``x_lim``, ``y_lim``, ``label_*``, ``repel_*``, colours,
    ...) go to the R function. ``metadata["genes"]`` has each gene's ``side``
    (``other``, ``sig_mid``, ``condA_ind`` / ``condA_rep`` / ``condB_ind`` /
    ``condB_rep``) and whether it is ``labelled``, once drawn.
    """
    missing = [c for c in RESPONSE_SPECIFICITY_COLUMNS if c not in combined]
    if missing:
        raise KeyError(f"Combined association table is missing column(s): {', '.join(missing)}.")
    if isinstance(conditions, dict):
        conditions = (conditions["condA"], conditions["condB"])
    if len(conditions) != 2:
        raise ValueError("conditions must name exactly two contexts.")
    data = combined.loc[:, list(RESPONSE_SPECIFICITY_COLUMNS)].apply(pd.to_numeric, errors="coerce")
    data.insert(0, "row_id", combined.index.astype(str))
    data = data.reset_index(drop=True)
    rendered, preview, value = _draw(
        "chronopack_plot_response_specificity_volcano_csv",
        output,
        tables={"data": data},
        condition_key={"condA": str(conditions[0]), "condB": str(conditions[1])},
        width=float(width),
        height=float(height),
        dpi=int(dpi),
        **kwargs,
    )
    genes = _r_to_pandas(value) if value is not None else None
    return RPlotResult(
        "response_specificity_volcano",
        data=data,
        output=rendered,
        preview_png=preview,
        metadata={"conditions": tuple(conditions), "n_genes": int(len(data)), "genes": genes},
    )

__all__ = [
    "CATEGORICAL_COLORS",
    "LOW_SIGNIFICANCE_COLOR",
    "MODE_COLORS",
    "MODE_ORDER",
    "RPlotResult",
    "association_volcano",
    "cluster_medoids_and_top_genes",
    "clustering_metrics",
    "fdr_wald_break",
    "gene_trajectory",
    "response_specificity_volcano",
    "trajectory_features",
]
