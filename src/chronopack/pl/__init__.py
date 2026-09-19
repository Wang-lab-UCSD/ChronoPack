"""Plotting wrappers for ChronoPack results."""

from __future__ import annotations

from chronopack._backends.r import get_r_resource_path, read_r_resource_text


def r_plotting_resource_path(name="dimred_plotting_helpers.R"):
    """Yield a path to a bundled R plotting helper resource."""
    return get_r_resource_path(name)


def read_r_plotting_resource(name="dimred_plotting_helpers.R"):
    """Read a bundled R plotting helper resource."""
    return read_r_resource_text(name)


def __getattr__(name):
    if name in {"RPlotResult", "gene_trajectory", "cluster_trajectories", "trajectory_features"}:
        from . import _gene_trajectory_plots

        return getattr(_gene_trajectory_plots, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "RPlotResult",
    "cluster_trajectories",
    "gene_trajectory",
    "r_plotting_resource_path",
    "read_r_plotting_resource",
    "trajectory_features",
]
