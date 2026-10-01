"""Plotting wrappers for gene-expression trajectory AnnData outputs."""

from __future__ import annotations

from dataclasses import dataclass
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


def _render(plot_function: str, data: pd.DataFrame, output, **kwargs) -> str | None:
    if output is None:
        return None
    output = str(Path(output))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("R-backed plotting requires optional dependency 'rpy2'.") from exc
    source_r_resource("dimred_plotting_helpers.R")
    with TemporaryDirectory() as tmpdir:
        data_path = Path(tmpdir) / "plot_data.csv"
        data.to_csv(data_path, index=False)
        r[plot_function](str(data_path), output, **kwargs)
    return output


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
    rendered = _render(
        "chronopack_plot_gene_trajectory_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
        normalize=bool(normalize),
    )
    return RPlotResult("gene_trajectory", data=data, output=rendered, metadata={"normalize": normalize})


def cluster_trajectories(
    adata,
    *,
    groupby: str = "broad_cluster",
    fit_key: str = GENE_TRAJECTORY_KEY,
    clustering_key: str = GENE_TRAJECTORY_CLUSTERING_KEY,
    rows=None,
    normalize: bool = True,
    output: str | Path | None = None,
    width: float = 8,
    height: float = 5,
    title: str | None = None,
) -> RPlotResult:
    """Plot median +/- IQR fitted trajectories for dynamic-mode groups."""
    del fit_key, clustering_key
    selected = pd.Index(adata.obs_names.astype(str)) if rows is None else pd.Index([str(x) for x in rows])
    long_df = _trajectory_long(adata, rows=selected, normalize=normalize)
    if groupby not in long_df.columns:
        raise KeyError(f"Grouping column '{groupby}' is not available in the trajectory plot data.")
    data = (
        long_df.groupby([groupby, "time"], dropna=False)
        .agg(
            median=("value", "median"),
            q25=("value", lambda x: x.quantile(0.25)),
            q75=("value", lambda x: x.quantile(0.75)),
            n=("row_id", "nunique"),
        )
        .reset_index()
        .rename(columns={groupby: "cluster"})
    )
    data["cluster"] = data["cluster"].astype(str)
    rendered = _render(
        "chronopack_plot_cluster_trajectories_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
        normalize=bool(normalize),
    )
    return RPlotResult("cluster_trajectories", data=data, output=rendered, metadata={"groupby": groupby})


def trajectory_features(
    adata,
    *,
    x: str = "PC1",
    y: str = "PC2",
    color: str = "broad_cluster",
    basis: str | None = None,
    clustering_key: str = GENE_TRAJECTORY_CLUSTERING_KEY,
    fit_key: str = GENE_TRAJECTORY_KEY,
    output: str | Path | None = None,
    width: float = 6,
    height: float = 5,
    title: str | None = None,
) -> RPlotResult:
    """Plot two trajectory features against each other, coloured by ``color``.

    ``x`` and ``y`` are ``traj.obs`` columns (descriptors), or with ``basis`` the
    components (``"PC1"``, ``"PC2"``, ...) of ``traj.obsm[basis]``, e.g. the
    ``X_pca`` coordinates stored by the clustering notebook; trajectories
    without coordinates there (NaN) are left out.
    """
    del clustering_key, fit_key
    data = _gene_labels(adata).reset_index(drop=True)
    for name in dict.fromkeys([x, y, color] if basis is None else [color]):
        if name not in adata.obs:
            raise KeyError(f"Trajectory observation column '{name}' is not available.")
        data[name] = adata.obs[name].to_numpy()
    if basis is not None:
        if basis not in adata.obsm:
            raise KeyError(f"Trajectory embedding adata.obsm['{basis}'] is not available.")
        coords = np.asarray(adata.obsm[basis], dtype=float)
        for name in (x, y):
            data[name] = coords[:, int(str(name).removeprefix("PC")) - 1]
        data = data.dropna(subset=[x, y]).reset_index(drop=True)
    rendered = _render(
        "chronopack_plot_trajectory_features_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
        x=str(x),
        y=str(y),
        color=str(color),
    )
    return RPlotResult("trajectory_features", data=data, output=rendered, metadata={"x": x, "y": y, "color": color})


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


__all__ = ["RPlotResult", "gene_trajectory", "cluster_trajectories", "trajectory_features", "clustering_metrics"]
