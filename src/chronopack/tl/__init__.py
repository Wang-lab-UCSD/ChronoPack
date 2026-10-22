"""Analysis tools for ChronoPack."""

__all__ = [
    "c2c",
    "grn",
    "CellStateTrajectoryResult",
    "GeneDynamicsResult",
    "GeneTrajectoryClusteringResult",
    "get_trajectories",
    "ot_force_directed",
    "cluster_cells_to_trajectories",
    "compute_terminal_fate_probabilities",
    "coupling_from_kernel",
    "fit_gene_dynamics",
    "predict_gene_dynamics",
    "observed_gene_dynamics",
    "prescreen_time_association",
    "signed_max_log2fc_from_start",
    "combine_context_associations",
    "early_fate_correlations",
    "rotate_embedding_to_time",
    "identify_tag_bias_candidates",
    "filter_tag_bias_candidates",
    "tag_bias_exclusion_genes",
    "merge_trajectory_graphs",
    "trajectory_soft_memberships",
    "align_gene_trajectories",
    "compute_gene_trajectory_descriptors",
    "cluster_gene_trajectories",
    "standardize_gene_trajectories",
    "topk_mode_recovery",
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
    if name in {
        "fit_gene_dynamics",
        "predict_gene_dynamics",
        "observed_gene_dynamics",
        "prescreen_time_association",
        "signed_max_log2fc_from_start",
        "combine_context_associations",
        "identify_tag_bias_candidates",
        "filter_tag_bias_candidates",
        "tag_bias_exclusion_genes",
    }:
        from . import _gene_dynamics

        return getattr(_gene_dynamics, name)
    if name in {
        "align_gene_trajectories",
            "compute_gene_trajectory_descriptors",
        "cluster_gene_trajectories",
            "standardize_gene_trajectories",
    }:
        from . import _gene_clustering

        return getattr(_gene_clustering, name)
    if name == "topk_mode_recovery":
        from . import _mode_conservation

        return getattr(_mode_conservation, name)
    if name == "get_trajectories":
        from ._trajectories import get_trajectories

        return get_trajectories
    if name in {
        "cluster_cells_to_trajectories",
        "cluster_cell_state_trajectories",
        "compute_terminal_fate_probabilities",
        "coupling_from_kernel",
        "early_fate_correlations",
        "merge_trajectory_graphs",
        "ot_force_directed",
        "rotate_embedding_to_time",
        "trajectory_soft_memberships",
    }:
        from . import _trajectories

        return getattr(_trajectories, name)
    if name in {"c2c", "grn"}:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
