"""Cell-cell communication analysis tools."""

from ._helpers import build_pair_series, split_by_cell_type
from ._inference import C2CResult, infer

__all__ = ["C2CResult", "build_pair_series", "infer", "split_by_cell_type"]
