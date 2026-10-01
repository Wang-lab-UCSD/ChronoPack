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
    if name in {
        "RPlotResult",
        "MODE_COLORS",
        "MODE_ORDER",
        "CATEGORICAL_COLORS",
        "LOW_SIGNIFICANCE_COLOR",
        "gene_trajectory",
        "trajectory_features",
        "cluster_medoids_and_top_genes",
        "clustering_metrics",
    }:
        from . import _gene_trajectory_plots

        return getattr(_gene_trajectory_plots, name)
    if name in {
        "mode_conservation",
        "mode_proportions",
        "mode_sankey",
        "mode_chord",
        "cluster_pair_representative_genes",
        "cluster_violin",
    }:
        from . import _mode_conservation_plots

        return getattr(_mode_conservation_plots, name)
    if name == "c2c":
        import importlib

        return importlib.import_module(f"{__name__}.c2c")
    if name == "cell_trajectories":
        import importlib

        return importlib.import_module(f"{__name__}.cell_trajectories")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "CATEGORICAL_COLORS",
    "LOW_SIGNIFICANCE_COLOR",
    "MODE_COLORS",
    "MODE_ORDER",
    "RPlotResult",
    "cluster_medoids_and_top_genes",
    "cluster_pair_representative_genes",
    "cluster_violin",
    "clustering_metrics",
    "gene_trajectory",
    "mode_chord",
    "mode_conservation",
    "mode_proportions",
    "mode_sankey",
    "c2c",
    "cell_trajectories",
    "r_plotting_resource_path",
    "read_r_plotting_resource",
    "trajectory_features",
]
