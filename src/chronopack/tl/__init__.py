"""Analysis tools for ChronoPack."""

__all__ = [
    "c2c",
    "grn",
    "CellStateTrajectoryResult",
    "GeneDynamicsResult",
    "GeneTrajectoryAlignmentResult",
    "GeneTrajectoryCharacterizationResult",
    "GeneTrajectoryClusteringResult",
    "GeneTrajectoryFinalizedClusteringResult",
    "ModeConservationResult",
    "get_trajectories",
    "cluster_cells_to_trajectories",
    "compute_terminal_fate_probabilities",
    "fit_gene_dynamics",
    "align_gene_trajectories",
    "compute_gene_trajectory_descriptors",
    "characterize_gene_trajectories",
    "cluster_gene_trajectories",
    "finalize_gene_trajectory_clusters",
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
    if name in {
        "GeneTrajectoryAlignmentResult",
        "GeneTrajectoryCharacterizationResult",
        "GeneTrajectoryClusteringResult",
        "GeneTrajectoryFinalizedClusteringResult",
    }:
        from . import _gene_clustering

        return getattr(_gene_clustering, name)
    if name == "ModeConservationResult":
        from ._mode_conservation import ModeConservationResult

        return ModeConservationResult
    if name == "fit_gene_dynamics":
        from ._gene_dynamics import fit_gene_dynamics

        return fit_gene_dynamics
    if name in {
        "align_gene_trajectories",
        "compute_gene_trajectory_descriptors",
        "characterize_gene_trajectories",
        "cluster_gene_trajectories",
        "finalize_gene_trajectory_clusters",
    }:
        from . import _gene_clustering

        return getattr(_gene_clustering, name)
    if name in {"mode_conservation", "mode_switchers"}:
        from . import _mode_conservation

        return getattr(_mode_conservation, name)
    if name == "get_trajectories":
        from ._trajectories import get_trajectories

        return get_trajectories
    if name in {"cluster_cells_to_trajectories", "cluster_cell_state_trajectories", "compute_terminal_fate_probabilities"}:
        from . import _trajectories

        return getattr(_trajectories, name)
    if name in {"c2c", "grn"}:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
