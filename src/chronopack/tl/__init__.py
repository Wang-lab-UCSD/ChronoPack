"""Analysis tools for ChronoPack."""

__all__ = [
    "c2c",
    "grn",
    "get_trajectories",
    "cluster_cell_state_trajectories",
    "compute_terminal_fate_probabilities",
    "fit_gene_dynamics",
    "characterize_gene_trajectories",
    "cluster_gene_trajectories",
]


def __getattr__(name):
    if name == "fit_gene_dynamics":
        from ._gene_dynamics import fit_gene_dynamics

        return fit_gene_dynamics
    if name in {"characterize_gene_trajectories", "cluster_gene_trajectories"}:
        from . import _gene_clustering

        return getattr(_gene_clustering, name)
    if name == "get_trajectories":
        from ._trajectories import get_trajectories

        return get_trajectories
    if name in {"cluster_cell_state_trajectories", "compute_terminal_fate_probabilities"}:
        from . import _trajectories

        return getattr(_trajectories, name)
    if name in {"c2c", "grn"}:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
