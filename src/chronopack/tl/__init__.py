"""Analysis tools for ChronoPack."""

__all__ = [
    "c2c",
    "grn",
    "CellStateTrajectoryResult",
    "GeneDynamicsResult",
    "GeneTrajectoryClusteringResult",
    "ModeConservationResult",
    "get_trajectories",
    "cluster_cells_to_trajectories",
    "compute_terminal_fate_probabilities",
    "fit_gene_dynamics",
    "identify_tag_bias_candidates",
    "filter_tag_bias_candidates",
    "tag_bias_exclusion_genes",
    "merge_trajectory_graphs",
    "align_gene_trajectories",
    "association_paths_from_dimred_roots",
    "compute_gene_trajectory_descriptors",
    "cluster_gene_trajectories",
    "sce_paths_from_dimred_roots",
    "standardize_gene_trajectories",
    "mode_conservation",
    "mode_switchers",
]


def __getattr__(name):
    if name == "CellStateTrajectoryResult":
        from ._trajectories import CellStateTrajectoryResult

        return CellStateTrajectoryResult
    if name == "GeneDynamicsResult":
        from ._gene_dynamics import GeneDynamicsResult

        return GeneDynamicsResult
    if name == "GeneTrajectoryClusteringResult":
        from . import _gene_clustering

        return getattr(_gene_clustering, name)
    if name == "ModeConservationResult":
        from ._mode_conservation import ModeConservationResult

        return ModeConservationResult
    if name in {"fit_gene_dynamics", "identify_tag_bias_candidates", "filter_tag_bias_candidates", "tag_bias_exclusion_genes"}:
        from . import _gene_dynamics

        return getattr(_gene_dynamics, name)
    if name in {
        "align_gene_trajectories",
        "association_paths_from_dimred_roots",
        "compute_gene_trajectory_descriptors",
        "cluster_gene_trajectories",
        "sce_paths_from_dimred_roots",
        "standardize_gene_trajectories",
    }:
        from . import _gene_clustering

        return getattr(_gene_clustering, name)
    if name in {"mode_conservation", "mode_switchers"}:
        from . import _mode_conservation

        return getattr(_mode_conservation, name)
    if name == "get_trajectories":
        from ._trajectories import get_trajectories

        return get_trajectories
    if name in {
        "cluster_cells_to_trajectories",
        "cluster_cell_state_trajectories",
        "compute_terminal_fate_probabilities",
        "merge_trajectory_graphs",
    }:
        from . import _trajectories

        return getattr(_trajectories, name)
    if name in {"c2c", "grn"}:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
