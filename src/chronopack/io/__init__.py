"""Input/output and AnnData preparation helpers."""

from __future__ import annotations

from chronopack._core import validate_adata, validate_preprocessed_adata
from ._preprocessing import (
    TIME_GROUP_CATEGORIES,
    add_time_metadata,
    annotate_qc_metrics,
    filter_cells_by_qc,
    filter_low_count_genes,
    filter_technical_genes,
    initialize_expression_layers,
    load_time_tagged_anndata,
    map_time_to_group,
    read_time_tagged_dge,
    score_cell_cycle,
)


def read_h5ad(path, **kwargs):
    """Read an AnnData object."""
    import anndata as ad

    return ad.read_h5ad(path, **kwargs)


def write_h5ad(adata, path, **kwargs):
    """Write an AnnData object after common ChronoPack validation."""
    validate_adata(adata)
    return adata.write_h5ad(path, **kwargs)


validate = validate_adata

__all__ = [
    "TIME_GROUP_CATEGORIES",
    "add_time_metadata",
    "annotate_qc_metrics",
    "filter_cells_by_qc",
    "filter_low_count_genes",
    "filter_technical_genes",
    "initialize_expression_layers",
    "load_time_tagged_anndata",
    "map_time_to_group",
    "read_time_tagged_dge",
    "read_h5ad",
    "score_cell_cycle",
    "validate",
    "validate_preprocessed_adata",
    "write_h5ad",
]


def write_tradeseq_inputs(*args, **kwargs):
    """Write tradeSeq inputs using the legacy-compatible helper.

    This helper is retained for migrated notebooks. New workflows should call
    :func:`chronopack.tl.fit_gene_dynamics`, which prepares tradeSeq inputs as
    an internal backend step.
    """
    from chronopack._legacy.tradeseq import WriteTradeSeqInputsFromAnnData

    return WriteTradeSeqInputsFromAnnData(*args, **kwargs)
