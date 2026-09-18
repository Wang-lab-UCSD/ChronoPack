"""Input/output and AnnData preparation helpers."""

from __future__ import annotations

from chronopack._core import validate_adata


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
    "annotate_qc_flags",
    "filter_cells_by_qc",
    "filter_low_count_genes",
    "filter_technical_genes",
    "load_time_tagged_anndata",
    "map_time_to_group",
    "read_h5ad",
    "score_cell_cycle",
    "validate",
    "write_h5ad",
    "write_tradeseq_inputs",
]


_PREPROCESSING_EXPORTS = {
    "TIME_GROUP_CATEGORIES",
    "add_time_metadata",
    "annotate_qc_flags",
    "filter_cells_by_qc",
    "filter_low_count_genes",
    "filter_technical_genes",
    "load_time_tagged_anndata",
    "map_time_to_group",
    "score_cell_cycle",
}


def write_tradeseq_inputs(*args, **kwargs):
    """Write tradeSeq inputs using the legacy-compatible helper."""
    from chronopack._legacy.tradeseq import WriteTradeSeqInputsFromAnnData

    return WriteTradeSeqInputsFromAnnData(*args, **kwargs)


def __getattr__(name):
    if name in _PREPROCESSING_EXPORTS:
        from chronopack._legacy import preprocessing

        return getattr(preprocessing, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
