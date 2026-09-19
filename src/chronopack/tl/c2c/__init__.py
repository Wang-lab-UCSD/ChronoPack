"""Cell-cell communication analysis tools."""

from chronopack._legacy.c2c import (  # noqa: F401
    all_expressed,
    build_pair_series,
    detection_series,
    filter_expressed_secreted_lr_catalog,
    load_secreted_cellchat_pairs,
    load_tradeseq_grid_exports,
    lr_activity,
    module_series,
    plot_avg_vs_logp,
    plot_genes_over_time,
    plot_pair_overlay,
    restrict_ligands_to_association,
    score_lr_pairs,
    split_by_cell_type,
    split_genes,
    time_series,
    time_to_minute,
    zscore,
)

from ._inference import C2CResult, infer

__all__ = ["C2CResult", "infer"]
