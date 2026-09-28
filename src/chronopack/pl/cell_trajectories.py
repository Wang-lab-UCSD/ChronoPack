"""Plotting helpers for cell-state trajectory and fate outputs."""

from __future__ import annotations

from chronopack._core._keys import CELL_STATE_TRAJECTORY_CLUSTER_KEY, OT_GRAPH_KEY, OT_UMAP_KEY
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


__all__ = ["fate_probabilities", "mass_retention", "proportions_over_time"]
