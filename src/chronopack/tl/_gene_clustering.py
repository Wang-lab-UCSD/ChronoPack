"""Gene-expression trajectory characterization and dynamic-mode clustering."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import (
    GENE_TRAJECTORY_CLUSTERING_KEY,
    GENE_TRAJECTORY_KEY,
)


@dataclass
class GeneTrajectoryAlignmentResult:
    """Result for aligned fitted trajectories and standardized matrices."""

    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    adata: Any | None = None
    pc_df: pd.DataFrame | None = None
    time_grid: np.ndarray | None = None
    sidecar_paths: dict[str, str] = field(default_factory=dict)
    backend_result: Any | None = None


@dataclass
class GeneTrajectoryCharacterizationResult:
    """Lightweight result for standardized fitted gene-expression trajectories."""

    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    adata: Any | None = None
    pc_df: pd.DataFrame | None = None
    time_grid: np.ndarray | None = None
    sidecar_paths: dict[str, str] = field(default_factory=dict)
    backend_result: Any | None = None


@dataclass
class GeneTrajectoryClusteringResult:
    """Lightweight result metadata for dynamic-mode clustering."""

    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    adata: Any | None = None
    pc_df: pd.DataFrame | None = None
    metrics: dict[str, pd.DataFrame] = field(default_factory=dict)
    selected_k: dict[str, dict[str, Any]] = field(default_factory=dict)
    sidecar_paths: dict[str, str] = field(default_factory=dict)
    backend_result: Any | None = None


@dataclass
class GeneTrajectoryFinalizedClusteringResult:
    """Combined induced/repressed cluster assignments."""

    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    adata: Any | None = None
    pc_df: pd.DataFrame | None = None
    sidecar_paths: dict[str, str] = field(default_factory=dict)
    backend_result: Any | None = None


def _read_matrix(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def _ensure_gene_adata(adata, var_names):
    import anndata as ad

    var_names = pd.Index(var_names.astype(str) if hasattr(var_names, "astype") else list(map(str, var_names)))
    if adata is None:
        return ad.AnnData(
            X=np.zeros((0, len(var_names)), dtype=np.float32),
            var=pd.DataFrame(index=var_names),
        )
    validate_adata(adata)
    if list(map(str, adata.var_names)) != var_names.tolist():
        current = pd.Index(map(str, adata.var_names))
        if current.is_unique and var_names.is_unique and set(current) == set(var_names):
            return adata[:, var_names].copy()
        raise ValueError("adata.var_names must match the gene-trajectory row order.")
    return adata


def _r_named_vector(values: dict[str, Any], *, numeric=True):
    from rpy2.robjects import FloatVector, StrVector

    keys = list(values.keys())
    vals = [values[k] for k in keys]
    vec = FloatVector([float(v) for v in vals]) if numeric else StrVector([str(v) for v in vals])
    vec.names = StrVector(keys)
    return vec


def _r_character_vector(values):
    from rpy2.robjects import StrVector

    return StrVector([str(x) for x in values])


def _r_int_vector(values):
    from rpy2.robjects import IntVector

    return IntVector([int(x) for x in values])


def _r_scalar_list(values: dict[str, Any]):
    from rpy2.robjects import BoolVector, FloatVector, IntVector, ListVector, StrVector

    out = {}
    for key, value in values.items():
        if isinstance(value, bool):
            out[str(key)] = BoolVector([value])
        elif isinstance(value, int) and not isinstance(value, bool):
            out[str(key)] = IntVector([value])
        elif isinstance(value, float):
            out[str(key)] = FloatVector([value])
        elif isinstance(value, str):
            out[str(key)] = StrVector([value])
        else:
            out[str(key)] = value
    return ListVector(out)


def _r_config(config):
    from rpy2.robjects import ListVector, StrVector

    out = {}
    for key, entry in config.items():
        out[str(key)] = ListVector(
            {
                "dir": str(entry["dir"]),
                "asso": StrVector([str(x) for x in entry["asso"]]),
                "sig": StrVector([str(x) for x in entry["sig"]]),
                "sce": StrVector([str(x) for x in entry["sce"]]),
            }
        )
    return ListVector(out)


def dimred_config(root_puro, root_k562a, root_lps):
    """Return the current fixed-sidecar configuration for DimensionalityReduction."""
    from rpy2.robjects import r

    source_r_resource("dimred_clustering_helpers.R")
    return r["chronopack_dimred_config"](str(root_puro), str(root_k562a), str(root_lps))


def _default_outputs(output_dir):
    output_dir = Path(output_dir)
    return {
        "fit": output_dir / "gene_trajectory_fit.csv",
        "time_grid": output_dir / "gene_trajectory_time_grid.csv",
        "pc_df_aligned": output_dir / "pc_df.aligned.csv",
        "pc_df_characterized": output_dir / "pc_df.characterized.csv",
        "fit_z": output_dir / "gene_trajectory_fit_z.csv",
        "fit_ddt_z": output_dir / "gene_trajectory_fit_ddt_z.csv",
        "fit_significant_z": output_dir / "gene_trajectory_fit_significant_z.csv",
        "fit_ddt_significant_z": output_dir / "gene_trajectory_fit_ddt_significant_z.csv",
        "pc_df": output_dir / "pc_df.new.csv",
        "pc_df_rds": output_dir / "pc_df.new.rds",
        "x_pca": output_dir / "gene_trajectory_X_pca.csv",
        "xd_pca": output_dir / "gene_trajectory_Xd_pca.csv",
        "induced_metrics": output_dir / "cluster_metrics.induced.csv",
        "repressed_metrics": output_dir / "cluster_metrics.repressed.csv",
    }


def _write_characterization_to_adata(adata, *, key_added, pc_df, fit, time_grid, ddt=None):
    for col in pc_df.columns:
        adata.var[f"{key_added}_{col}"] = pc_df[col].to_numpy()
    adata.varm[f"{key_added}_fit"] = fit.loc[adata.var_names].to_numpy()
    if ddt is not None:
        adata.varm[f"{key_added}_fit_ddt"] = ddt.loc[adata.var_names].to_numpy()
    record_provenance(
        adata,
        key_added,
        analysis="gene_expression_trajectory_characterization",
        parameters={"time_grid": [float(x) for x in time_grid]},
        result_keys={
            "var": [f"{key_added}_{col}" for col in pc_df.columns],
            "varm": [f"{key_added}_fit"] + ([f"{key_added}_fit_ddt"] if ddt is not None else []),
        },
    )


def _write_clustering_to_adata(adata, *, key_added, pc_df, x_pca=None):
    for col in pc_df.columns:
        if col not in adata.var:
            adata.var[f"{key_added}_{col}"] = pc_df[col].to_numpy()
    if x_pca is not None:
        adata.varm[f"{key_added}_X_pca"] = x_pca.loc[adata.var_names].to_numpy()
    cluster_cols = [
        col
        for col in [
            "cluster_broad_induced",
            "cluster_broad_repressed",
            "cluster_fine_induced",
            "cluster_fine_repressed",
            "broad_cluster",
            "cluster_dist_broad",
            "cluster_dist_fine",
        ]
        if col in pc_df.columns
    ]
    record_provenance(
        adata,
        key_added,
        analysis="gene_expression_trajectory_clustering",
        parameters={},
        result_keys={
            "var": [f"{key_added}_{col}" for col in cluster_cols if f"{key_added}_{col}" in adata.var],
            "varm": [f"{key_added}_X_pca"] if x_pca is not None else [],
        },
    )


def align_gene_trajectories(
    adata=None,
    *,
    config=None,
    root_puro=None,
    root_k562a=None,
    root_lps=None,
    output_dir,
    key_added=GENE_TRAJECTORY_KEY,
    keys=("K562a", "K562_puro", "K562b", "THP1"),
    t_interval_map=None,
    time_map=None,
    offset_map=None,
    scale_map=None,
    pvalue_threshold=0.05,
    lineage_pick=1,
    copy=False,
    return_model=True,
):
    """Align fitted trajectories across experiments and compute standardized matrices."""
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("Gene-trajectory alignment requires optional dependency 'rpy2'.") from exc

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_r_resource("dimred_clustering_helpers.R")

    if config is None:
        if root_puro is None or root_k562a is None or root_lps is None:
            raise ValueError("Provide either config or root_puro/root_k562a/root_lps.")
        r_config = r["chronopack_dimred_config"](str(root_puro), str(root_k562a), str(root_lps))
    else:
        r_config = _r_config(config)

    t_interval_map = t_interval_map or {"K562a": 10, "K562_puro": 10, "K562b": 15, "THP1": 15}
    time_map = time_map or {"K562a": 110, "K562_puro": 110, "K562b": 165, "THP1": 165}
    offset_map = offset_map or {"K562a": 0, "K562_puro": 0, "K562b": -68, "THP1": -10}
    scale_map = scale_map or {"THP1": 0.65}

    backend_result = r["chronopack_align_gene_trajectories"](
        r_config,
        str(output_dir),
        keys=_r_character_vector(keys),
        t_interval_map=_r_named_vector(t_interval_map),
        time_map=_r_named_vector(time_map),
        offset_map=_r_named_vector(offset_map),
        scale_map=_r_named_vector(scale_map),
        pvalue_threshold=float(pvalue_threshold),
        lineage_pick=int(lineage_pick),
        write_outputs=True,
    )

    paths = _default_outputs(output_dir)
    fit = _read_matrix(paths["fit"])
    pc_df = pd.read_csv(paths["pc_df_aligned"], index_col=0)
    time_grid = pd.read_csv(paths["time_grid"])["time"].to_numpy(dtype=float)
    ddt = _read_matrix(paths["fit_ddt_z"])

    target = _ensure_gene_adata(adata.copy() if copy and adata is not None else adata, pc_df.index)
    _write_characterization_to_adata(target, key_added=key_added, pc_df=pc_df, fit=fit, time_grid=time_grid, ddt=ddt)

    sidecars = {name: str(path) for name, path in paths.items() if path.exists()}
    metadata = {
        "output_dir": str(output_dir),
        "keys": list(keys),
        "t_interval_map": dict(t_interval_map),
        "time_map": dict(time_map),
        "offset_map": dict(offset_map),
        "scale_map": dict(scale_map),
        "pvalue_threshold": pvalue_threshold,
        "lineage_pick": lineage_pick,
        "time_grid": [float(x) for x in time_grid],
    }
    result = GeneTrajectoryAlignmentResult(
        key_added=key_added,
        metadata=metadata,
        adata=target,
        pc_df=pc_df,
        time_grid=time_grid,
        sidecar_paths=sidecars,
        backend_result=backend_result,
    )
    if copy:
        return (target, result) if return_model else target
    return result if return_model else None


def compute_gene_trajectory_descriptors(
    adata=None,
    *,
    alignment,
    output_dir,
    key_added=GENE_TRAJECTORY_KEY,
    pvalue_threshold=None,
    descriptor_args=None,
    copy=False,
    return_model=True,
):
    """Compute trajectory descriptors from an aligned trajectory result."""
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("Gene-trajectory descriptor calculation requires optional dependency 'rpy2'.") from exc

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_r_resource("dimred_clustering_helpers.R")

    kwargs = {
        "alignment": alignment.backend_result,
        "output_dir": str(output_dir),
        "write_outputs": True,
    }
    if pvalue_threshold is not None:
        kwargs["pvalue_threshold"] = float(pvalue_threshold)
    if descriptor_args is not None:
        kwargs["descriptor_args"] = _r_scalar_list(descriptor_args)

    backend_result = r["chronopack_compute_gene_trajectory_descriptors"](**kwargs)

    paths = _default_outputs(output_dir)
    fit = _read_matrix(paths["fit"])
    pc_df = pd.read_csv(paths["pc_df_characterized"], index_col=0)
    time_grid = pd.read_csv(paths["time_grid"])["time"].to_numpy(dtype=float)
    ddt = _read_matrix(paths["fit_ddt_z"])

    base_adata = adata if adata is not None else alignment.adata
    target = _ensure_gene_adata(base_adata.copy() if copy and base_adata is not None else base_adata, pc_df.index)
    _write_characterization_to_adata(target, key_added=key_added, pc_df=pc_df, fit=fit, time_grid=time_grid, ddt=ddt)

    sidecars = {name: str(path) for name, path in paths.items() if path.exists()}
    metadata = dict(alignment.metadata)
    metadata.update(
        {
            "output_dir": str(output_dir),
            "pvalue_threshold": pvalue_threshold if pvalue_threshold is not None else alignment.metadata.get("pvalue_threshold"),
            "time_grid": [float(x) for x in time_grid],
        }
    )
    result = GeneTrajectoryCharacterizationResult(
        key_added=key_added,
        metadata=metadata,
        adata=target,
        pc_df=pc_df,
        time_grid=time_grid,
        sidecar_paths=sidecars,
        backend_result=backend_result,
    )
    if copy:
        return (target, result) if return_model else target
    return result if return_model else None


def characterize_gene_trajectories(
    adata=None,
    *,
    config=None,
    root_puro=None,
    root_k562a=None,
    root_lps=None,
    output_dir,
    key_added=GENE_TRAJECTORY_KEY,
    keys=("K562a", "K562_puro", "K562b", "THP1"),
    t_interval_map=None,
    time_map=None,
    offset_map=None,
    scale_map=None,
    pvalue_threshold=0.05,
    lineage_pick=1,
    copy=False,
    return_model=True,
):
    """Compatibility wrapper: align trajectories, then compute descriptors."""
    alignment = align_gene_trajectories(
        adata=adata,
        config=config,
        root_puro=root_puro,
        root_k562a=root_k562a,
        root_lps=root_lps,
        output_dir=output_dir,
        key_added=key_added,
        keys=keys,
        t_interval_map=t_interval_map,
        time_map=time_map,
        offset_map=offset_map,
        scale_map=scale_map,
        pvalue_threshold=pvalue_threshold,
        lineage_pick=lineage_pick,
        copy=copy,
        return_model=True,
    )
    return compute_gene_trajectory_descriptors(
        adata=alignment.adata,
        alignment=alignment,
        output_dir=output_dir,
        key_added=key_added,
        pvalue_threshold=pvalue_threshold,
        copy=copy,
        return_model=return_model,
    )


def cluster_gene_trajectories(
    adata=None,
    *,
    characterization=None,
    characterization_output_dir=None,
    output_dir,
    rows,
    cluster_name,
    key_added=GENE_TRAJECTORY_CLUSTERING_KEY,
    alpha=None,
    k_broad=None,
    k_fine=None,
    k_grid=range(2, 101),
    wss_fraction=0.95,
    variance_threshold=0.99,
    pca_columns=4,
    copy=False,
    return_model=True,
):
    """Cluster an explicitly supplied subset of standardized gene trajectories."""
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("Gene-trajectory clustering requires optional dependency 'rpy2'.") from exc

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_r_resource("dimred_clustering_helpers.R")
    if alpha is None:
        raise ValueError("Pass alpha explicitly.")

    if characterization is None:
        if characterization_output_dir is None:
            raise ValueError("Provide characterization result or characterization_output_dir.")
        # Reconstruct the small R characterization object from standardized sidecars.
        # This path is intended for notebook restarts and regression scripts.
        raise ValueError("Restarting clustering from sidecar-only characterization is not implemented yet.")

    backend_kwargs = {
        "row_ids": _r_character_vector(rows),
        "cluster_name": str(cluster_name),
        "alpha": float(alpha),
        "K_grid": _r_int_vector(k_grid),
        "wss_fraction": float(wss_fraction),
        "variance_threshold": float(variance_threshold),
        "pca_columns": int(pca_columns),
        "write_outputs": True,
    }
    if k_broad is not None or k_fine is not None:
        if k_broad is None or k_fine is None:
            raise ValueError("Pass both k_broad and k_fine, or omit both to evaluate metrics only.")
        backend_kwargs["k_broad"] = int(k_broad)
        backend_kwargs["k_fine"] = int(k_fine)

    backend_result = r["chronopack_cluster_gene_trajectories"](
        characterization.backend_result,
        str(output_dir),
        **backend_kwargs,
    )

    paths = _default_outputs(output_dir)
    cluster_pc_df_path = output_dir / f"pc_df.{cluster_name}.csv"
    pc_df = pd.read_csv(cluster_pc_df_path, index_col=0)
    x_pca = _read_matrix(paths["x_pca"]) if paths["x_pca"].exists() else None
    metrics = {}
    metrics_path = output_dir / f"cluster_metrics.{cluster_name}.csv"
    if metrics_path.exists():
        metrics[str(cluster_name)] = pd.read_csv(metrics_path)

    target = _ensure_gene_adata(adata.copy() if copy and adata is not None else adata, pc_df.index)
    _write_clustering_to_adata(target, key_added=key_added, pc_df=pc_df, x_pca=x_pca)

    sidecars = {name: str(path) for name, path in paths.items() if path.exists()}
    metadata = {
        "output_dir": str(output_dir),
        "cluster_name": str(cluster_name),
        "alpha": alpha,
        "rows": list(map(str, rows)),
        "k_broad": None if k_broad is None else int(k_broad),
        "k_fine": None if k_fine is None else int(k_fine),
        "k_grid": list(map(int, k_grid)),
        "wss_fraction": wss_fraction,
        "variance_threshold": variance_threshold,
        "pca_columns": pca_columns,
    }
    result = GeneTrajectoryClusteringResult(
        key_added=key_added,
        metadata=metadata,
        adata=target,
        pc_df=pc_df,
        metrics=metrics,
        sidecar_paths=sidecars,
        backend_result=backend_result,
    )
    if copy:
        return (target, result) if return_model else target
    return result if return_model else None


def finalize_gene_trajectory_clusters(
    adata=None,
    *,
    characterization,
    clusterings,
    output_dir,
    key_added=GENE_TRAJECTORY_CLUSTERING_KEY,
    cell_type_order=("K562a", "K562b", "K562_puro", "THP1"),
    copy=False,
    return_model=True,
):
    """Combine explicitly computed trajectory cluster assignments."""
    try:
        from rpy2.robjects import ListVector, r
    except ImportError as exc:
        raise ImportError("Gene-trajectory clustering requires optional dependency 'rpy2'.") from exc

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_r_resource("dimred_clustering_helpers.R")

    r_clusterings = ListVector({res.metadata["cluster_name"]: res.backend_result for res in clusterings})
    backend_result = r["chronopack_finalize_gene_trajectory_clusters"](
        characterization.backend_result,
        r_clusterings,
        str(output_dir),
        cell_type_order=_r_character_vector(cell_type_order),
        write_outputs=True,
    )

    paths = _default_outputs(output_dir)
    pc_df = pd.read_csv(paths["pc_df"], index_col=0)
    x_pca = _read_matrix(paths["x_pca"]) if paths["x_pca"].exists() else None
    base_adata = adata if adata is not None else characterization.adata
    target = _ensure_gene_adata(base_adata.copy() if copy and base_adata is not None else base_adata, pc_df.index)
    _write_clustering_to_adata(target, key_added=key_added, pc_df=pc_df, x_pca=x_pca)

    sidecars = {name: str(path) for name, path in paths.items() if path.exists()}
    metadata = {
        "output_dir": str(output_dir),
        "cluster_names": [res.metadata["cluster_name"] for res in clusterings],
        "cell_type_order": list(cell_type_order),
    }
    result = GeneTrajectoryFinalizedClusteringResult(
        key_added=key_added,
        metadata=metadata,
        adata=target,
        pc_df=pc_df,
        sidecar_paths=sidecars,
        backend_result=backend_result,
    )
    if copy:
        return (target, result) if return_model else target
    return result if return_model else None


__all__ = [
    "GeneTrajectoryAlignmentResult",
    "GeneTrajectoryCharacterizationResult",
    "GeneTrajectoryClusteringResult",
    "GeneTrajectoryFinalizedClusteringResult",
    "align_gene_trajectories",
    "compute_gene_trajectory_descriptors",
    "characterize_gene_trajectories",
    "cluster_gene_trajectories",
    "finalize_gene_trajectory_clusters",
    "dimred_config",
]
