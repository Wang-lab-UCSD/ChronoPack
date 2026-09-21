"""Compatibility imports for migrated gene-dynamics helpers."""

from chronopack._utils._deprecation import warn_legacy_module
from chronopack.tl._gene_dynamics import (  # noqa: F401
    compute_smooth_correlations,
    filter_tag_bias_candidates,
    identify_tag_bias_candidates,
    load_negative_control_logcpm,
    prepare_tradeseq_inputs,
    refilter_tag_bias_candidates,
    smooth_series,
    tag_bias_exclusion_genes,
    temporal_pseudobulk_from_anndata,
)

warn_legacy_module("chronopack.tradeseq", "chronopack.tl.fit_gene_dynamics")

TemporalPseudobulkFromAnnData = temporal_pseudobulk_from_anndata
ComputeSmoothCorrelations = compute_smooth_correlations


def WriteTradeSeqInputsFromAnnData(tmp, key=None, tag=""):
    """Compatibility wrapper preserving the historical notebook signature."""
    if key is None:
        raise ValueError("key must be provided")
    return prepare_tradeseq_inputs(tmp, f"../intmd/{key}", tag=tag)

__all__ = [
    "ComputeSmoothCorrelations",
    "TemporalPseudobulkFromAnnData",
    "WriteTradeSeqInputsFromAnnData",
    "compute_smooth_correlations",
    "filter_tag_bias_candidates",
    "identify_tag_bias_candidates",
    "load_negative_control_logcpm",
    "prepare_tradeseq_inputs",
    "refilter_tag_bias_candidates",
    "smooth_series",
    "tag_bias_exclusion_genes",
    "temporal_pseudobulk_from_anndata",
]
