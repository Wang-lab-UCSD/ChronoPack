"""AnnData schema conventions for current ChronoPack workflows."""

from __future__ import annotations

from dataclasses import dataclass

from ._keys import DEFAULT_COUNTS_LAYER, DEFAULT_TIME_KEY

TIME_TAG_KEY = "FINAL TIME TAG"
TIME_GROUP_KEY = "time_group"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
DATASET_KEY = "dataset"
CELL_TYPE_KEY = "cell_type"

CPM_LAYER = "cpm"
LOG1P_LAYER = "log1p"
CORRECTED_LAYER = "corrected"

MITO_FLAG_KEY = "mt"
RIBOSOMAL_FLAG_KEY = "rb"

TIME_GROUP_CATEGORIES = (
    "Group 1 (T1-3)",
    "Group 2 (T4-6)",
    "Group 3 (T7-9)",
    "Group 4 (T10-12)",
)


@dataclass(frozen=True)
class PreprocessingSchema:
    """Minimal AnnData contract consumed by current downstream workflows."""

    time_key: str = DEFAULT_TIME_KEY
    time_group_key: str = TIME_GROUP_KEY
    sample_key: str = SAMPLE_KEY
    counts_layer: str = DEFAULT_COUNTS_LAYER
    cpm_layer: str = CPM_LAYER
    log1p_layer: str = LOG1P_LAYER
    corrected_layer: str = CORRECTED_LAYER
    cell_type_key: str = CELL_TYPE_KEY
    mito_flag_key: str = MITO_FLAG_KEY
    ribosomal_flag_key: str = RIBOSOMAL_FLAG_KEY


PREPROCESSING_SCHEMA = PreprocessingSchema()


def map_time_to_group(time_tag) -> str:
    """Map integer time tags 1..12 onto the existing four coarse groups."""
    time_tag = int(time_tag)
    if 1 <= time_tag <= 3:
        return TIME_GROUP_CATEGORIES[0]
    if 4 <= time_tag <= 6:
        return TIME_GROUP_CATEGORIES[1]
    if 7 <= time_tag <= 9:
        return TIME_GROUP_CATEGORIES[2]
    if 10 <= time_tag <= 12:
        return TIME_GROUP_CATEGORIES[3]
    return "Unknown Group"
