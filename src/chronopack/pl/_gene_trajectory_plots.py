"""Plotting wrappers for standardized gene-expression trajectory outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core._keys import (
    CHRONOPACK_UNS_KEY,
    GENE_TRAJECTORY_CLUSTERING_KEY,
    GENE_TRAJECTORY_KEY,
)


@dataclass
class RPlotResult:
    """Small result object for R-backed ChronoPack plots."""

    plot_type: str
    data: pd.DataFrame
    output: str | None = None
    backend: str = "R/ggplot2"
    metadata: dict[str, Any] | None = None


def _prefixed_var_col(adata, key: str, col: str) -> str:
    candidates = [f"{key}_{col}", col]
    for candidate in candidates:
        if candidate in adata.var.columns:
            return candidate
    raise KeyError(f"Could not find AnnData var column '{key}_{col}' or '{col}'.")


def _fit_matrix(adata, fit_key: str) -> np.ndarray:
    key = f"{fit_key}_fit"
    if key not in adata.varm:
        raise KeyError(f"Missing standardized fitted trajectory matrix adata.varm['{key}'].")
    fit = np.asarray(adata.varm[key], dtype=float)
    if fit.shape[0] != adata.n_vars:
        raise ValueError(f"adata.varm['{key}'] has {fit.shape[0]} rows but adata has {adata.n_vars} genes.")
    return fit


def _time_grid(adata, fit_key: str, n_timepoints: int) -> np.ndarray:
    entry = adata.uns.get(CHRONOPACK_UNS_KEY, {}).get(fit_key, {})
    params = entry.get("parameters", {}) if isinstance(entry, dict) else {}
    values = params.get("time_grid")
    if values is None:
        return np.arange(n_timepoints, dtype=float)
    out = np.asarray(values, dtype=float)
    if out.shape[0] != n_timepoints:
        raise ValueError(
            f"Stored time grid for '{fit_key}' has length {out.shape[0]}, expected {n_timepoints}."
        )
    return out


def _gene_labels(adata, fit_key: str) -> pd.DataFrame:
    gene_col = _prefixed_var_col(adata, fit_key, "gene") if f"{fit_key}_gene" in adata.var else None
    cell_type_col = _prefixed_var_col(adata, fit_key, "cell_type") if f"{fit_key}_cell_type" in adata.var else None
    labels = pd.DataFrame(index=adata.var_names.astype(str))
    labels["row_id"] = labels.index
    labels["gene"] = adata.var[gene_col].astype(str).to_numpy() if gene_col else labels["row_id"].str.replace(r"_(?=[^_]+$).*", "", regex=True)
    labels["cell_type"] = adata.var[cell_type_col].astype(str).to_numpy() if cell_type_col else ""
    return labels


def _select_rows(adata, *, fit_key: str, genes=None) -> pd.Index:
    var_names = pd.Index(adata.var_names.astype(str))
    if genes is None:
        return var_names
    requested = pd.Index([str(g) for g in genes])
    direct = requested[requested.isin(var_names)]
    remaining = requested.difference(direct)
    selected = list(direct)
    if len(remaining):
        labels = _gene_labels(adata, fit_key)
        gene_matches = labels.index[labels["gene"].isin(remaining)].tolist()
        selected.extend(gene_matches)
    missing = requested.difference(pd.Index(selected).union(pd.Index(labels["gene"] if "labels" in locals() else [])))
    if len(selected) == 0 or len(missing) > 0:
        raise KeyError(f"Could not find requested genes/rows: {', '.join(missing.astype(str))}")
    return pd.Index(selected).drop_duplicates()


def _z_row(matrix: np.ndarray) -> np.ndarray:
    mu = np.nanmean(matrix, axis=1, keepdims=True)
    sd = np.nanstd(matrix, axis=1, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return (matrix - mu) / sd


def _trajectory_long(
    adata,
    *,
    fit_key: str,
    clustering_key: str | None,
    rows: pd.Index,
    normalize: bool,
) -> pd.DataFrame:
    fit = _fit_matrix(adata, fit_key)
    row_lookup = pd.Series(np.arange(adata.n_vars), index=adata.var_names.astype(str))
    idx = row_lookup.loc[rows].to_numpy()
    values = fit[idx, :]
    if normalize:
        values = _z_row(values)
    time = _time_grid(adata, fit_key, values.shape[1])

    labels = _gene_labels(adata, fit_key).loc[rows].reset_index(drop=True)
    meta = labels.copy()
    for col in ["waldStat_1", "peak_log2fc", "t_onset", "t_peak"]:
        try:
            meta[col] = adata.var.loc[rows, _prefixed_var_col(adata, fit_key, col)].to_numpy()
        except KeyError:
            pass
    if clustering_key is not None:
        for col in ["broad_cluster", "cluster_fine_induced", "cluster_fine_repressed"]:
            try:
                meta[col] = adata.var.loc[rows, _prefixed_var_col(adata, clustering_key, col)].astype(str).to_numpy()
            except KeyError:
                pass

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
    """Plot standardized fitted trajectories for selected genes or gene-context rows."""
    rows = _select_rows(adata, fit_key=fit_key, genes=genes)
    data = _trajectory_long(adata, fit_key=fit_key, clustering_key=clustering_key, rows=rows, normalize=normalize)
    rendered = _render(
        "chronopack_plot_gene_trajectory_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
        normalize=bool(normalize),
    )
    return RPlotResult("gene_trajectory", data=data, output=rendered, metadata={"fit_key": fit_key})


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
    """Plot median +/- IQR fitted trajectories for gene dynamic-mode groups."""
    selected = pd.Index(adata.var_names.astype(str)) if rows is None else pd.Index([str(x) for x in rows])
    long_df = _trajectory_long(adata, fit_key=fit_key, clustering_key=clustering_key, rows=selected, normalize=normalize)
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
    return RPlotResult(
        "cluster_trajectories",
        data=data,
        output=rendered,
        metadata={"fit_key": fit_key, "clustering_key": clustering_key, "groupby": groupby},
    )


def trajectory_features(
    adata,
    *,
    x: str = "PC1",
    y: str = "PC2",
    color: str = "broad_cluster",
    clustering_key: str = GENE_TRAJECTORY_CLUSTERING_KEY,
    fit_key: str = GENE_TRAJECTORY_KEY,
    output: str | Path | None = None,
    width: float = 6,
    height: float = 5,
    title: str | None = None,
) -> RPlotResult:
    """Plot gene trajectory features, descriptors, or PCA coordinates."""
    data = _gene_labels(adata, fit_key).reset_index(drop=True)

    def add_feature(name: str) -> None:
        pca_key = f"{clustering_key}_X_pca"
        if name.startswith("PC") and pca_key in adata.varm:
            pos = int(name[2:]) - 1
            data[name] = np.asarray(adata.varm[pca_key])[:, pos]
            return
        try:
            col = _prefixed_var_col(adata, clustering_key, name)
        except KeyError:
            col = _prefixed_var_col(adata, fit_key, name)
        data[name] = adata.var[col].to_numpy()

    for name in dict.fromkeys([x, y, color]):
        add_feature(name)

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
    return RPlotResult(
        "trajectory_features",
        data=data,
        output=rendered,
        metadata={"fit_key": fit_key, "clustering_key": clustering_key, "x": x, "y": y, "color": color},
    )


__all__ = ["RPlotResult", "gene_trajectory", "cluster_trajectories", "trajectory_features"]
