"""Plotting helpers for cell-state trajectory and fate outputs."""

from __future__ import annotations

from chronopack._core._keys import CELL_STATE_TRAJECTORY_CLUSTER_KEY, OT_GRAPH_KEY, OT_UMAP_KEY
from chronopack.tl._trajectories import (
    PlotHardClusterMassProportionsOverTime,
    mass_retention_plots_and_umap,
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


__all__ = ["mass_retention", "proportions_over_time"]
