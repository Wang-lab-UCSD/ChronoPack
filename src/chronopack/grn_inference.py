"""Compatibility imports for ChronoPack GRN inference helpers."""

from chronopack._utils._deprecation import warn_legacy_module

warn_legacy_module("chronopack.grn_inference", "chronopack.tl.grn")

from chronopack._legacy.grn_inference import load_reference_tfs, load_valid_targets
from chronopack.tl.grn import (  # noqa: F401
    GRNResult,
    TargetGeneScoreResult,
    build_chronopack_ranked_edges,
    infer,
    rank_regulators,
    score_targets,
    write_chronopack_ranked_edges,
)

__all__ = [
    "GRNResult",
    "TargetGeneScoreResult",
    "infer",
    "rank_regulators",
    "score_targets",
    "build_chronopack_ranked_edges",
    "write_chronopack_ranked_edges",
    "load_valid_targets",
    "load_reference_tfs",
]
