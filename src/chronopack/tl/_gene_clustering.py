"""Gene-expression trajectory characterization and dynamic-mode clustering."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from chronopack._core import record_provenance, validate_adata


@dataclass
class GeneTrajectoryClusteringResult:
    """Lightweight result metadata for fitted gene-expression curve clustering."""

    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)


def characterize_gene_trajectories(
    adata,
    *,
    fit_key="cp_fit",
    key_added="gene_trajectory_features",
    copy=False,
    **kwargs,
):
    """Record characterization of fitted gene-expression trajectories."""
    target = adata.copy() if copy else adata
    validate_adata(target)
    record_provenance(
        target,
        key_added,
        analysis="gene_expression_trajectory_features",
        parameters={"fit_key": fit_key, **kwargs},
    )
    return target if copy else None


def cluster_gene_trajectories(
    adata,
    *,
    fit_key="cp_fit",
    key_added="gene_dynamic_modes",
    copy=False,
    **kwargs,
):
    """Record clustering of fitted gene-expression trajectories into modes."""
    target = adata.copy() if copy else adata
    validate_adata(target)
    record_provenance(
        target,
        key_added,
        analysis="gene_expression_trajectory_clustering",
        parameters={"fit_key": fit_key, **kwargs},
    )
    result = GeneTrajectoryClusteringResult(key_added=key_added, metadata={"fit_key": fit_key, **kwargs})
    if copy:
        return target, result
    return result


__all__ = [
    "GeneTrajectoryClusteringResult",
    "characterize_gene_trajectories",
    "cluster_gene_trajectories",
]
