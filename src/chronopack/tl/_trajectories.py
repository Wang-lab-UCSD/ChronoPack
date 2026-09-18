"""Cell-state trajectory and fate analysis public wrappers."""

from __future__ import annotations

from chronopack._core import record_provenance, validate_adata


def get_trajectories(
    adata,
    *,
    tmap_dir,
    time_key="time",
    key_added="cell_state_trajectories",
    copy=False,
    **kwargs,
):
    """Compute WOT-derived cell-state trajectory outputs for an AnnData object."""
    from chronopack._legacy.trajectories import run_wot_transport_workflow

    target = adata.copy() if copy else adata
    validate_adata(target, time_key=time_key)
    result = run_wot_transport_workflow(target, tmap_dir=tmap_dir, **kwargs)
    record_provenance(
        target,
        key_added,
        analysis="cell_state_trajectories",
        backend="waddington_ot",
        parameters={"time_key": time_key, **kwargs},
        result_keys={
            "obsp": ["ot_graph"],
            "obsm": ["X_ot_umap"],
            "sidecar": ["tmaps", "transitions_long", "tmats"],
        },
    )
    return (target, result) if copy else result


__all__ = [
    "CalculateTerminalFatesAcrossCells",
    "build_transition_matrices",
    "cosine_topk_leiden_from_terminal_dist",
    "extract_wot_couplings",
    "get_trajectories",
    "read_tmap",
    "run_wot_transport_workflow",
    "sparse_adjacency_from_tmats",
    "terminal_fate_trajectory_clustering",
    "terminal_fates_from_terminal_clusters",
    "topk_csr_rows",
]


_LEGACY_TRAJECTORY_EXPORTS = set(__all__) - {"get_trajectories"}


def __getattr__(name):
    if name in _LEGACY_TRAJECTORY_EXPORTS:
        from chronopack._legacy import trajectories

        return getattr(trajectories, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
