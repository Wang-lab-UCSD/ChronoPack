"""Input/output and AnnData preparation helpers."""

from __future__ import annotations

from chronopack._core import validate_adata
from ._arrayexpress import download_arrayexpress_files
from ._preprocessing import (
    TIME_GROUP_CATEGORIES,
    add_time_metadata,
    filter_low_count_genes,
    filter_technical_genes,
    load_time_tagged_anndata,
    map_time_to_group,
    read_time_tagged_dge,
    score_cell_cycle,
    sum_dge_by_time_tag,
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
    "download_arrayexpress_files",
    "TIME_GROUP_CATEGORIES",
    "add_time_metadata",
    "filter_low_count_genes",
    "filter_technical_genes",
    "load_time_tagged_anndata",
    "map_time_to_group",
    "read_time_tagged_dge",
    "read_h5ad",
    "score_cell_cycle",
    "validate",
    "write_h5ad",
    "sum_dge_by_time_tag",
]
