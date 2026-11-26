"""Gene-expression trajectory AnnData construction and clustering."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import (
    CHRONOPACK_UNS_KEY,
    GENE_TRAJECTORY_CLUSTERING_KEY,
    GENE_TRAJECTORY_KEY,
)


@dataclass
class GeneTrajectoryClusteringResult:
    """Candidate clustering solutions and diagnostics for one trajectory subset."""

    metrics: pd.DataFrame
    assignments: dict[int, pd.Series]
    centroid_distances: dict[int, pd.Series] = field(default_factory=dict)
    pca: pd.DataFrame | None = None
    # Features (rows = trajectories) the distances are computed on: the joint shape + derivative features, or
    # their top ``n_pcs`` principal-component scores when ``n_pcs`` is set.
    features: pd.DataFrame | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
    key_added: str = GENE_TRAJECTORY_CLUSTERING_KEY


def _read_matrix(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


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


def _r_sce_dict(sce_dict: dict[str, Any]):
    from rpy2.robjects import ListVector, StrVector

    out = {}
    for key, value in sce_dict.items():
        if isinstance(value, (str, Path)):
            out[str(key)] = StrVector([str(value)])
        else:
            out[str(key)] = value
    return ListVector(out)


def _r_optional_path_dict(values: dict[str, Any] | None):
    if values is None:
        return None
    from rpy2.robjects import ListVector, StrVector

    return ListVector({str(k): StrVector([str(v)]) for k, v in values.items() if v is not None})


def _make_traj_adata(
    *,
    fit: pd.DataFrame,
    obs: pd.DataFrame,
    time_grid: np.ndarray,
    raw_payload: dict[str, dict[str, Any]],
    raw_size_bytes: int,
    raw_embedded: bool,
    metadata: dict[str, Any],
    key_added: str,
):
    import anndata as ad

    obs = obs.loc[fit.index].copy()
    obs.index = obs.index.astype(str)
    var = pd.DataFrame({"time": time_grid}, index=[f"t{i}" for i in range(fit.shape[1])])
    traj = ad.AnnData(X=fit.to_numpy(dtype=float), obs=obs, var=var)
    if CHRONOPACK_UNS_KEY not in traj.uns:
        traj.uns[CHRONOPACK_UNS_KEY] = {}
    traj.uns[CHRONOPACK_UNS_KEY]["alignment"] = {
        **metadata,
        "raw_trajectory_size_bytes": int(raw_size_bytes),
        "raw_trajectories_embedded": bool(raw_embedded),
    }
    if raw_embedded:
        traj.uns[CHRONOPACK_UNS_KEY]["raw_trajectories"] = raw_payload
    record_provenance(
        traj,
        key_added,
        analysis="gene_expression_trajectory_alignment",
        backend="R/tradeSeq",
        parameters=metadata,
        result_keys={"X": "aligned fitted expression", "obs": list(obs.columns), "var": ["time"]},
    )
    return traj


def align_gene_trajectories(
    sce_dict: dict[str, Any],
    *,
    key_added: str = GENE_TRAJECTORY_KEY,
    keys=None,
    time_map=None,
    offset_map=None,
    scale_map=None,
    association_paths: dict[str, Any] | None = None,
    lineage_pick: int = 1,
    embed_raw_max_bytes: int = 100_000_000,
):
    """Build a gene-context x aligned-time AnnData object from fitted tradeSeq SCEs."""
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("Gene-trajectory alignment requires optional dependency 'rpy2'.") from exc

    time_map = time_map or {"K562a": 110, "K562b": 165, "THP1": 165}
    offset_map = offset_map or {"K562a": 0, "K562b": -68, "THP1": -10}
    scale_map = scale_map or {"THP1": 0.65}
    # Contexts: every sce_dict entry, in order, unless a subset is requested.
    keys = tuple(sce_dict) if keys is None else tuple(k for k in keys if k in sce_dict)
    if not keys:
        raise ValueError("No requested keys are present in sce_dict.")

    source_r_resource("dimred_clustering_helpers.R")
    tmp_ctx = TemporaryDirectory()  # R hand-off files only
    work_dir = Path(tmp_ctx.name)
    try:
        # rpy2 cannot convert None; omit association_paths so R uses its NULL default.
        optional = {}
        if association_paths is not None:
            optional["association_paths"] = _r_optional_path_dict(association_paths)
        r["chronopack_align_sce_trajectories_to_files"](
            _r_sce_dict({k: sce_dict[k] for k in keys}),
            str(work_dir),
            keys=_r_character_vector(keys),
            time_map=_r_named_vector({k: time_map[k] for k in keys}),
            offset_map=_r_named_vector({k: offset_map[k] for k in keys}),
            scale_map=_r_named_vector(scale_map),
            lineage_pick=int(lineage_pick),
            write_raw=True,
            **optional,
        )
        fit = _read_matrix(work_dir / "trajectory_aligned_fit.csv")
        obs = pd.read_csv(work_dir / "trajectory_obs.csv", index_col=0)
        time_grid = pd.read_csv(work_dir / "trajectory_time_grid.csv")["time"].to_numpy(dtype=float)

        raw_payload: dict[str, dict[str, Any]] = {}
        raw_size_bytes = 0
        raw_meta_path = work_dir / "raw_fitted_trajectories.meta.csv"
        if raw_meta_path.exists():
            raw_meta = pd.read_csv(raw_meta_path)
            raw_size_bytes = int(raw_meta["size_bytes"].sum()) if "size_bytes" in raw_meta else 0
        raw_embedded = raw_size_bytes <= int(embed_raw_max_bytes)
        if raw_embedded:
            raw_dir = work_dir / "raw_fitted_trajectories"
            for key in keys:
                raw_fit = _read_matrix(raw_dir / f"{key}.fit.csv")
                raw_time = pd.read_csv(raw_dir / f"{key}.time.csv")["time"].to_numpy(dtype=float)
                raw_payload[key] = {
                    "fit": raw_fit.to_numpy(dtype=float),
                    "time": raw_time,
                    "obs_names": raw_fit.index.astype(str).to_numpy(),
                    "gene": raw_fit.index.astype(str).to_numpy(),
                    "context": key,
                }

        metadata = {
            "keys": list(keys),
            "time_map": {str(k): float(time_map[k]) for k in keys},
            "offset_map": {str(k): float(offset_map[k]) for k in keys},
            "scale_map": {str(k): float(v) for k, v in dict(scale_map).items()},
            "scale_map_use": {
                str(k): float(dict(scale_map).get(k, 1.0))
                for k in keys
            },
            "source_sce": {str(k): str(sce_dict[k]) for k in keys if isinstance(sce_dict[k], (str, Path))},
            "association_paths_used_for_migration": {
                str(k): str(v) for k, v in dict(association_paths or {}).items()
            },
            "lineage_pick": int(lineage_pick),
            "aligned_time_grid": [float(x) for x in time_grid],
        }
        traj = _make_traj_adata(
            fit=fit,
            obs=obs,
            time_grid=time_grid,
            raw_payload=raw_payload,
            raw_size_bytes=raw_size_bytes,
            raw_embedded=raw_embedded,
            metadata=metadata,
            key_added=key_added,
        )
    finally:
        tmp_ctx.cleanup()
    return traj


def compute_gene_trajectory_descriptors(
    traj,
    *,
    inplace: bool = True,
    descriptor_args: dict[str, Any] | None = None,
    key_added: str = GENE_TRAJECTORY_KEY,
):
    """Compute trajectory descriptors from ``traj.X`` and append them to ``traj.obs``."""
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("Gene-trajectory descriptor calculation requires optional dependency 'rpy2'.") from exc

    validate_adata(traj)
    target = traj if inplace else traj.copy()
    time_grid = np.asarray(target.var["time"], dtype=float)
    source_r_resource("dimred_clustering_helpers.R")

    tmp_ctx = TemporaryDirectory()  # R hand-off files only
    work_dir = Path(tmp_ctx.name)
    try:
        traj_path = work_dir / "descriptor_input_fit.csv"
        time_path = work_dir / "descriptor_input_time.csv"
        desc_path = work_dir / "trajectory_descriptors.csv"
        pd.DataFrame(target.X, index=target.obs_names, columns=target.var_names).to_csv(traj_path)
        pd.DataFrame({"time": time_grid}).to_csv(time_path, index=False)
        alignment_meta = target.uns.get(CHRONOPACK_UNS_KEY, {}).get("alignment", {})
        scale_map = alignment_meta.get("scale_map_use", alignment_meta.get("scale_map", {}))
        if descriptor_args is None:
            descriptor_args = {
                "N_minutes": 10 / len(time_grid),
                "z_diff_threshold": 0.2,
                "require_consecutive": 30,
                "onset_mode": "horizontal",
                "pseudocount": 0,
            }
        r_kwargs = {}
        if scale_map:
            r_kwargs["scale_map_use"] = _r_named_vector(scale_map)
        if descriptor_args:
            r_kwargs["descriptor_args"] = _r_scalar_list(descriptor_args)
        r["chronopack_compute_descriptors_from_files"](
            str(traj_path),
            str(time_path),
            str(desc_path),
            **r_kwargs,
        )
        desc = pd.read_csv(desc_path, index_col=0)
        desc = desc.loc[target.obs_names]
    finally:
        tmp_ctx.cleanup()

    for col in desc.columns:
        target.obs[col] = desc[col].to_numpy()
    record_provenance(
        target,
        f"{key_added}_descriptors",
        analysis="gene_expression_trajectory_descriptors",
        backend="R",
        parameters={"descriptor_args": descriptor_args or {}},
        result_keys={"obs": list(desc.columns)},
    )
    return None if inplace else target


def _write_subset_matrix(traj, work_dir: Path, basename: str) -> tuple[Path, Path]:
    fit_path = work_dir / f"{basename}.fit.csv"
    time_path = work_dir / f"{basename}.time.csv"
    pd.DataFrame(traj.X, index=traj.obs_names, columns=traj.var_names).to_csv(fit_path)
    pd.DataFrame({"time": np.asarray(traj.var["time"], dtype=float)}).to_csv(time_path, index=False)
    return fit_path, time_path


def cluster_gene_trajectories(
    traj,
    *,
    k=range(3, 20),
    alpha: float,
    reference=None,
    key_added: str = GENE_TRAJECTORY_CLUSTERING_KEY,
    wss_fraction: float = 0.95,
    variance_threshold: float = 0.99,
    pca_columns: int = 4,
    n_pcs: int | None = None,
):
    """Compute candidate clustering labels and metrics for every requested k once.

    Trajectories are clustered by Ward.D2 hierarchical clustering on Euclidean distances between joint
    features: row-z-scored fitted curves and row-z-scored time derivatives, each block scaled by its total
    column variance and weighted by ``sqrt(1 - alpha)`` and ``sqrt(alpha)``.

    Parameters
    ----------
    traj
        Gene-trajectory AnnData from :func:`align_gene_trajectories` (obs = trajectories, var = time points).
    k
        Candidate numbers of clusters; one ``cutree`` of the same tree per value.
    alpha
        Weight of the derivative block, in [0, 1].
    reference
        Optional larger trajectory AnnData (same time points, containing every row of ``traj``). Standardization,
        block scaling and the ``n_pcs`` PCA are fit on it and ``traj`` is placed within it.
    wss_fraction
        WSS-resolution fraction used by the k-selection metrics.
    variance_threshold, pca_columns
        Settings of the separate descriptive PCA of the curves (``result.pca``). They do not affect the
        clustering.
    n_pcs
        ``None`` (default): cluster on the full joint features. An integer: cluster on the top ``n_pcs``
        principal-component scores of the joint features (centred, not scaled, not whitened; fit on
        ``reference`` when given). The k-selection metrics and centroid distances are computed in the same
        space. ``n_pcs`` equal to the rank of the joint features reproduces the full-feature clustering.
    """
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("Gene-trajectory clustering requires optional dependency 'rpy2'.") from exc

    validate_adata(traj)
    if alpha is None:
        raise ValueError("Pass alpha explicitly.")
    k_grid = [int(x) for x in k]
    if not k_grid:
        raise ValueError("At least one candidate k is required.")
    if n_pcs is not None:
        if isinstance(n_pcs, bool) or int(n_pcs) != n_pcs or int(n_pcs) < 1:
            raise ValueError("n_pcs must be None or a positive integer.")
        n_pcs = int(n_pcs)

    source_r_resource("dimred_clustering_helpers.R")
    tmp_ctx = TemporaryDirectory()  # R hand-off files only
    work_dir = Path(tmp_ctx.name)
    try:
        fit_path, time_path = _write_subset_matrix(traj, work_dir, "cluster_input")
        reference_path = None
        if reference is not None:
            if list(map(str, reference.var_names)) != list(map(str, traj.var_names)):
                raise ValueError("reference.var_names must match traj.var_names for clustering.")
            reference_path, _ = _write_subset_matrix(reference, work_dir, "cluster_reference")
        r["chronopack_cluster_trajectory_matrix_to_files"](
            str(fit_path),
            str(time_path),
            str(work_dir),
            alpha=float(alpha),
            K_grid=_r_int_vector(k_grid),
            wss_fraction=float(wss_fraction),
            variance_threshold=float(variance_threshold),
            pca_columns=int(pca_columns),
            **({"reference_traj_mat_path": str(reference_path)} if reference_path is not None else {}),
            **({"n_pcs": n_pcs} if n_pcs is not None else {}),
        )
        metrics = pd.read_csv(work_dir / "cluster_metrics.csv")
        if "K" in metrics.columns and "k" not in metrics.columns:
            metrics.insert(0, "k", metrics["K"].astype(int))
        assignments_df = _read_matrix(work_dir / "cluster_assignments.csv")
        distances_df = _read_matrix(work_dir / "cluster_centroid_distances.csv")
        pca = _read_matrix(work_dir / "gene_trajectory_X_pca.csv")
        features = _read_matrix(work_dir / "cluster_features.csv").loc[traj.obs_names.astype(str)]
        feature_pca = pd.read_csv(work_dir / "cluster_feature_pca.csv")
    finally:
        tmp_ctx.cleanup()

    expected_index = pd.Index(traj.obs_names.astype(str))
    assignments: dict[int, pd.Series] = {}
    centroid_distances: dict[int, pd.Series] = {}
    for kk in k_grid:
        col = f"k{kk}"
        labels = pd.Series(assignments_df.loc[expected_index, col].astype(int).to_numpy(), index=expected_index, name=col)
        dists = pd.Series(distances_df.loc[expected_index, col].astype(float).to_numpy(), index=expected_index, name=col)
        if not labels.index.equals(expected_index):
            raise RuntimeError(f"Assignment index for k={kk} does not match traj.obs_names.")
        assignments[kk] = labels
        centroid_distances[kk] = dists

    metadata = {
        "alpha": float(alpha),
        "k_grid": k_grid,
        "wss_fraction": float(wss_fraction),
        "variance_threshold": float(variance_threshold),
        "pca_columns": int(pca_columns),
        "n_pcs": n_pcs,
        "n_pcs_variance_explained": (
            None if n_pcs is None else float(feature_pca["variance_explained"].iloc[0])
        ),
        "n_trajectories": int(traj.n_obs),
        "reference_n_trajectories": None if reference is None else int(reference.n_obs),
    }
    return GeneTrajectoryClusteringResult(
        metrics=metrics,
        assignments=assignments,
        centroid_distances=centroid_distances,
        pca=pca,
        features=features,
        metadata=metadata,
        key_added=key_added,
    )


def standardize_gene_trajectories(traj):
    """Return fitted, row-standardized and derivative-standardized trajectories.

    ``traj`` is the gene-trajectory AnnData from :func:`align_gene_trajectories`
    (obs = gene-context trajectories, var = aligned time points). Returns
    ``(fit, fit_z, ddt_z)`` DataFrames: the fitted curves, their row z-scores
    (sample SD; constant rows keep SD 1), and z-scores of the forward-difference
    time derivative, whose columns are the time points after the first.
    """

    def _z_row(M):
        M = M.astype(float)
        mu = M.mean(axis=1, keepdims=True)
        sd = M.std(axis=1, ddof=1, keepdims=True)
        sd[(~pd.notna(sd)) | (sd == 0)] = 1
        return (M - mu) / sd

    X = np.asarray(traj.X)
    time = traj.var["time"].to_numpy(dtype=float)
    ddt = (X[:, 1:] - X[:, :-1]) / (time[1:] - time[:-1])
    fit = pd.DataFrame(X, index=traj.obs_names, columns=traj.var_names)
    fit_z = pd.DataFrame(_z_row(X.copy()), index=traj.obs_names, columns=traj.var_names)
    ddt_z = pd.DataFrame(_z_row(ddt), index=traj.obs_names, columns=traj.var_names[1:])
    return fit, fit_z, ddt_z


__all__ = [
    "GeneTrajectoryClusteringResult",
    "align_gene_trajectories",
    "cluster_gene_trajectories",
    "compute_gene_trajectory_descriptors",
    "standardize_gene_trajectories",
]
